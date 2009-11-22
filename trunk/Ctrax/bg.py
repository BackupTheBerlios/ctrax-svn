# bg.py
# KMB 09/11/07

import numpy as num
from numpy import fft
from numpy import random
from scipy.linalg.basic import eps
import scipy.io
import sys
import threading
import time
import wx
from wx import xrc
from colormapk import colormap_image
import scipy.ndimage.measurements as meas
import scipy.ndimage.morphology as morph
from ellipsesk import find_ellipses
from ellipsesk import draw_ellipses
import imagesk
from params import params
import setarena
import fixbg
import roi

import motmot.wxvideo.wxvideo as wxvideo
#import motmot.wxglvideo.simple_overlay as wxvideo # part of Motmot
#import wxvalidatedtext as wxvt # part of Motmot
import motmot.wxvalidatedtext.wxvalidatedtext as wxvt # part of Motmot

import highboostfilter
#import histgmm
#from histgmm import make_vect
#from settings import const as settings_const
from version import DEBUG

import pkg_resources # part of setuptools
THRESH_RSRC_FILE = pkg_resources.resource_filename( __name__, "bg_thresh.xrc" )
SETTINGS_RSRC_FILE = pkg_resources.resource_filename( __name__, "bg_settings.xrc" )
HF_RSRC_FILE = pkg_resources.resource_filename( __name__, "homomorphic.xrc" )

SAVE_STUFF = False
USE_SAVED = False # USE_SAVED takes precedence over SAVE_STUFF

class NoMoreFramesException( Exception ):
    pass

class HomoFilt:
    """Implements a homomorphic filter and uses it to retrieve filtered
    movie frames."""
    def __init__( self , movie_size ):
        self.movie_size = movie_size
        self.hbf = highboostfilter.highboostfilter( self.movie_size,
                                                    params.hm_cutoff,
                                                    params.hm_order,
                                                    params.hm_boost )
        self.frame = None

    def SetFrame(self,parent):
        rsrc = xrc.XmlResource( HF_RSRC_FILE )
        self.frame = rsrc.LoadFrame( parent, "homomorphic_settings_frame" )
        self.hm_cutoff_box = xrc.XRCCTRL( self.frame, "hm_cutoff_input" )
        self.hm_boost_box = xrc.XRCCTRL( self.frame, "hm_boost_input" )
        self.hm_order_box = xrc.XRCCTRL( self.frame, "hm_order_input" )

        self.hm_cutoff_box.SetValue( str(params.hm_cutoff) )
        self.hm_boost_box.SetValue( str(params.hm_boost) )
        self.hm_order_box.SetValue( str(params.hm_order) )

        wxvt.setup_validated_float_callback( self.hm_cutoff_box,
                                             xrc.XRCID("hm_cutoff_input"),
                                             self.OnTextValidated,
                                             pending_color=params.wxvt_bg )
        wxvt.setup_validated_integer_callback( self.hm_boost_box,
                                               xrc.XRCID("hm_boost_input"),
                                               self.OnTextValidated,
                                               pending_color=params.wxvt_bg )
        wxvt.setup_validated_integer_callback( self.hm_order_box,
                                               xrc.XRCID("hm_order_input"),
                                               self.OnTextValidated,
                                               pending_color=params.wxvt_bg )

    def OnTextValidated( self, evt ):
        """Set global constants based on input settings."""
        # get new values
        new_cutoff = float( self.hm_cutoff_box.GetValue() )
        new_boost = float( self.hm_boost_box.GetValue() )
        new_order = int( self.hm_order_box.GetValue() )
        # validate ranges
        if new_cutoff >= 0.1 and new_cutoff <= 0.5 and \
               new_boost >= 1. and new_boost <= 10. and \
               new_order >= 2 and new_order <= 10:
            # update
            params.hm_cutoff = new_cutoff
            params.hm_boost = new_boost
            params.hm_order = new_order
            self.changed = True
        else:
            # revert
            self.hm_cutoff_box.SetValue( str(params.hm_cutoff) )
            self.hm_boost_box.SetValue( str(params.hm_boost) )
            self.hm_order_box.SetValue( str(params.hm_order) )

    def apply(self,I):
        Ifloat = I.astype( num.float )
        # compute fourier transform
        FFTlogI = num.fft.fft2( num.log( Ifloat + 1. ) )
        # apply high-boost filter
        img = num.exp( num.real( num.fft.ifft2( FFTlogI * self.hbf ) ) )
        # make sure it is in [0,1]
        img[img > 255.] = 255.
        img[img < 0.] = 0.
        return img

class BackgroundCalculator:
    def __init__( self, movie,
                  hm_cutoff=None, hm_boost=None, hm_order=None, use_thresh=None, use_median=False ):
        """Initialize background subtractor."""
        params.movie = movie
        #params = ConstantsBG() # so constants are attached if self is passed around

        params.n_frames = params.movie.get_n_frames()
        
        params.bg_firstframe = max(min(params.bg_firstframe,params.n_frames-1),0)
        params.bg_lastframe = max(params.bg_firstframe,params.bg_lastframe)

        # for now, hard code in everything
        #n_bg_frames = 10

        #params.n_bg_frames = n_bg_frames

        params.movie_size = (params.movie.get_height(), params.movie.get_width())
        params.GRID.setsize(params.movie_size)
        params.npixels = params.movie_size[0] * params.movie_size[1]

        self.show_img_type = params.SHOW_THRESH
        params.bg_type = params.BG_TYPE_OTHER
        params.bg_norm_type = params.BG_NORM_STD

        self.hf = HomoFilt(params.movie_size)

        # image of which parts flies are not able to walk in
        self.isarena = num.ones(params.movie_size,num.bool)
        self.roi = num.ones(params.movie_size,num.bool)

        if params.movie.type == 'sbfmf':
            self.center = params.movie.h_mov.bgcenter.copy()
            self.center.shape = params.movie.h_mov.framesize
            self.dev = params.movie.h_mov.bgstd.copy()
            self.dev.shape = params.movie.h_mov.framesize
            self.med = params.movie.h_mov.bgcenter.copy()
            self.med.shape = params.movie.h_mov.framesize
            self.mad = params.movie.h_mov.bgstd.copy()
            self.mad.shape = params.movie.h_mov.framesize
            self.mean = params.movie.h_mov.bgcenter.copy()
            self.mean.shape = params.movie.h_mov.framesize
            self.std = params.movie.h_mov.bgstd.copy()
            self.std.shape = params.movie.h_mov.framesize
            
            # approximate homomorphic filtering
            tmp = self.center.copy()
            issmall = tmp<1.
            tmp[issmall] = 1.
            self.hfnorm = self.hf.apply(self.center) / tmp
            self.hfnorm[issmall & (self.hfnorm<1.)] = 1.        

        # morphology stuff
        #print "creating morphology structures, with radii %d, %d"%(params.opening_radius,params.closing_radius)
        self.opening_struct = self.create_morph_struct(params.opening_radius)
        self.closing_struct = self.create_morph_struct(params.closing_radius)


    def updateParameters(self):

        if params.use_median:
            self.center = self.med.copy()
        else:
            self.center = self.mean.copy()

        if params.bg_norm_type == params.BG_NORM_STD:
            if params.use_median:
                self.dev = self.mad.copy()
            else:
                self.dev = self.std.copy()
            self.dev[self.dev < params.bg_std_min] = params.bg_std_min
            self.dev[self.dev > params.bg_std_max] = params.bg_std_max
        elif params.bg_norm_type == params.BG_NORM_INTENSITY:
            self.dev = self.center.copy()
        elif params.bg_norm_type == params.BG_NORM_HOMOMORPHIC:
            self.dev = self.hfnorm.copy()
        
        self.hf = HomoFilt(params.movie_size)

        # image of which parts flies are able to walk in, set using roi
        self.roi = roi.point_inside_polygons(params.movie_size[0],params.movie_size[1],
                                             params.roipolygons)

        self.UpdateIsArena()

    def meanstd( self, parent=None ):

        # shortcut names for oft-used values
        bg_lastframe = min(params.bg_lastframe,params.n_frames-1)
        nframes = bg_lastframe - params.bg_firstframe + 1
        n_bg_frames = min(nframes,params.n_bg_frames)

        if params.interactive and not params.batch_executing:

            # be prepared to cancel
            isbackup = False
            if hasattr(self,'mean'):
                mean0 = self.mean.copy()
                std0 = self.std.copy()
                isbackup = True

            progressbar = \
                wx.ProgressDialog('Computing Background Model',
                                  'Computing mean, standard deviation of %d frames to estimate background model'%n_bg_frames,
                                  n_bg_frames,
                                  parent,
                                  wx.PD_APP_MODAL|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT|wx.PD_REMAINING_TIME)

        nr = params.movie_size[0]
        nc = params.movie_size[1]
        # we will estimate background from every nframesskip-th frame
        nframesskip = int(num.floor(nframes/n_bg_frames))

        # initialize mean and std to 0
        self.mean = num.zeros((nr,nc))
        self.std = num.zeros((nr,nc))

        # main computation
        nframesused = 0
        for i in range(params.bg_firstframe,bg_lastframe+1,nframesskip):

            if params.interactive and not params.batch_executing:
                (keepgoing,skip) = progressbar.Update(value=nframesused+1,newmsg='Reading in frame %d (%d / %d)'%(i+1,nframesused,n_bg_frames))
                if not keepgoing:
                    progressbar.Destroy()

                    if isbackup:
                        self.mean = mean0
                        self.std = std0
                    else:
                        delattr(self,'mean')
                        delattr(self,'std')
                    
                    return False

            # read in the data
            im, stamp = params.movie.get_frame( int(i) )
            im = im.astype( num.float )
            # add to the mean
            self.mean += im
            # add to the variance
            self.std += im**2
            nframesused += 1

        # normalize
        self.mean /= num.double(nframesused)
        self.std /= num.double(nframesused)
        # actually compute variance, std
        self.std = num.sqrt(self.std - self.mean**2)
        
        if params.interactive and not params.batch_executing:
            progressbar.Destroy()

        return True

    def flexmedmad( self, parent=None ):
        bg_lastframe = min(params.bg_lastframe,params.n_frames-1)
        nframesmovie = bg_lastframe - params.bg_firstframe + 1

        nframes = min(params.n_bg_frames,nframesmovie)
        nr = params.movie_size[0]
        nc = params.movie_size[1]
        framesize = nr*nc
        # we will estimate background from every nframesskip-th frame
        #nframesskip = int(num.floor(10000/params.n_bg_frames))
        nframesskip = int(num.floor(nframesmovie/nframes))

        # which frame is the middle frame for computing the median?
        iseven = num.mod(nframes,2) == 0
        middle1 = num.floor(nframes/2)
        middle2 = middle1-1

        # number of rows to read in at a time; based on the assumption 
        # that we comfortably hold 100*(400x400) frames in memory. 
        nrsmall = num.int(num.floor(100.0*400.0*400.0/num.double(nc)/num.double(nframes)))
        if nrsmall < 1:
            nrsmall = 1
        if nrsmall > nr:
            nrsmall = nr
        # number of rows left in the last iteration might be less than nrsmall 
        nrsmalllast = num.mod(nr, nrsmall)
        # if evenly divides,  set last number of rows to nrsmall
        if nrsmalllast == 0:
            nrsmalllast = nrsmall
        # buffer holds the pixels for each frame that are read in. 
        # it holds nrsmall*nc pixels for each frame of nframes 
        buffersize = nrsmall*nc*nframes

        # prepare for cancel
        if params.interactive and not params.batch_executing:
            isbackup = False
            if hasattr(self,'med'):
                med0 = self.med.copy()
                mad0 = self.mad.copy()
                isbackup = True

            noffsets = num.ceil(float(nr) / float(nrsmall))

            progressbar = \
                wx.ProgressDialog('Computing Background Model',
                                  'Computing median, median absolute deviation of %d frames to estimate background model'%nframes,
                                  nframes*noffsets,
                                  parent,
                                  wx.PD_APP_MODAL|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT|wx.PD_REMAINING_TIME)


        # allocate memory for median and mad
        self.med = num.zeros((nr,nc))
        self.mad = num.zeros((nr,nc))

        offseti = 0

        for rowoffset in range(0,nr,nrsmall):
            if DEBUG: print "row offset = %d."%rowoffset

            # if this is the last iteration, there may not be npixelssmall left
            # store in npixelscurr the number of pixels that are to be read in,
            # store in seekperframecurr the amount to seek after reading.
            rowoffsetnext = int(min(rowoffset+nrsmall,nr))
            nrowscurr = rowoffsetnext - rowoffset

            # allocate memory for buffer that holds data read in in current pass
            buf = num.zeros((nframes,nrowscurr,nc),dtype=num.uint8)

            if DEBUG: print 'Reading ...'

            # loop through frames
            frame = params.bg_firstframe
            for i in range(nframes):

                if params.interactive and not params.batch_executing:
                    (keepgoing,skip) = progressbar.Update(value=offseti*nframes + i,newmsg='Reading in piece of frame %d (%d / %d), offset = %d (%d / %d)'%(frame,i+1,nframes,rowoffset,offseti+1,noffsets))
                    if not keepgoing:
                        progressbar.Destroy()
                        if isbackup:
                            self.med = med0
                            self.mad = mad0
                        else:
                            delattr(self,'med')
                            delattr(self,'mad')
                        return False

                # read in the entire frame
                data,stamp = params.movie.get_frame(frame)

                # crop out the rows we are interested in
                buf[i,:,:]  = data[rowoffset:rowoffsetnext,:]

                frame = frame + nframesskip

            # compute the median and median absolute difference at each 
            # pixel location
            if DEBUG: print 'Computing ...'
            # sort all the histories to get the median
            buf = num.rollaxis(num.rollaxis(buf,2,0),2,0)
            buf.sort(axis=2,kind='mergesort')
            # store the median
            self.med[rowoffset:rowoffsetnext,:] = buf[:,:,middle1]
            if iseven:
                self.med[rowoffset:rowoffsetnext,:] += buf[:,:,middle2]
                self.med[rowoffset:rowoffsetnext,:] /= 2.
                
            # store the absolute difference
            buf = num.double(buf)
            for j in range(nframes):
                buf[:,:,j] = num.abs(buf[:,:,j] - self.med[rowoffset:rowoffsetnext,:])

            # sort
            buf.sort(axis=2,kind='mergesort')

            # store the median absolute difference
            self.mad[rowoffset:rowoffsetnext,:] = buf[:,:,middle1]
            if iseven:
                self.mad[rowoffset:rowoffsetnext,:] += buf[:,:,middle2]
                self.mad[rowoffset:rowoffsetnext,:] /= 2.

            offseti += 1

        # estimate standard deviation assuming a Gaussian distribution
        # from the fact that half the data falls within mad
        # MADTOSTDFACTOR = 1./norminv(.75)
        MADTOSTDFACTOR = 1.482602
        self.mad *= MADTOSTDFACTOR

        if params.interactive and not params.batch_executing:
            progressbar.Destroy()
        return True
        
    def medmad( self, parent=None ):

        # shortcut names for oft-used values
        fp = params.movie.h_mov.file
        nr = params.movie_size[0]
        nc = params.movie_size[1]
 
        bg_lastframe = min(params.bg_lastframe,params.n_frames-1)
        nframesmovie = bg_lastframe - params.bg_firstframe + 1
        nframes = min(params.n_bg_frames,nframesmovie)

        if DEBUG: print "firstframe = %d, lastframe = %d"%(params.bg_firstframe,bg_lastframe)

       # we will estimate background from every nframesskip-th frame
        #nframesskip = int(num.floor(10000/params.n_bg_frames))
        nframesskip = int(num.floor(nframesmovie/nframes))

        bytesperchunk = params.movie.h_mov.bytes_per_chunk
        if num.mod(params.movie.h_mov.bits_per_pixel,8) != 0:
            print "Not sure what will happen if bytesperpixel is non-integer!!"
        bytesperpixel = params.movie.h_mov.bits_per_pixel/8.
        headersize = params.movie.h_mov.chunk_start+\
            params.bg_firstframe*bytesperchunk+\
            params.movie.h_mov.timestamp_len
        framesize = nr*nc
        nbytes = num.int(nr*nc*bytesperpixel)

        # which frame is the middle frame for computing the median?
        iseven = num.mod(nframes,2) == 0
        middle1 = num.floor(nframes/2)
        middle2 = middle1-1
        # which frame is the PROBWITHINONESTD th frame for computing mad
        #PROBWITHINONESTD = 0.6826894921371
        #madorder = num.double(nframes)*PROBWITHINONESTD
        #madorder1 = num.floor(madorder)
        #madorder2 = num.ceil(madorder)
        #if madorder1 == madorder2:
        #    madweight1 = .5
        #    madweight2 = .5
        #else:
        #    madweight1 = num.double(madorder2)-madorder
        #    madweight2 = madorder-num.double(madorder1)

        # number of rows to read in at a time; based on the assumption 
        # that we comfortably hold 100*(400x400) frames in memory. 
        nrsmall = num.int(num.floor(100.0*400.0*400.0/num.double(nc)/num.double(nframes)))
        if nrsmall < 1:
            nrsmall = 1
        if nrsmall > nr:
            nrsmall = nr
        # number of rows left in the last iteration might be less than nrsmall 
        nrsmalllast = num.mod(nr, nrsmall)
        # if evenly divides,  set last number of rows to nrsmall
        if nrsmalllast == 0:
            nrsmalllast = nrsmall
        # number of pixels corresponding to nrsmall rows
        nbytessmall = num.int(nrsmall*nc*bytesperpixel)
        # number of pixels corresponding to nrsmalllast rows
        nbytessmalllast = num.int(nrsmalllast*nc*bytesperpixel)
        # buffer holds the pixels for each frame that are read in. 
        # it holds npixelssmall for each frame of nframes 
        buffersize = nbytessmall*nframes
        # after we read in npixelssmall for a frame, we must seek forward
        # seekperframe 
        seekperframe = bytesperchunk*nframesskip-nbytessmall
        # in the last iteration, we only read in npixelssmalllast, so we
        # need to seek forward more
        seekperframelast = bytesperchunk*nframesskip-nbytessmalllast

        if params.interactive and not params.batch_executing:
            # be prepared to cancel
            isbackup = False
            if hasattr(self,'med'):
                med0 = self.med.copy()
                mad0 = self.mad.copy()
                isbackup = True

            noffsets = num.ceil(float(nbytes) / float(nbytessmall))

            progressbar = \
                wx.ProgressDialog('Computing Background Model',
                                  'Computing median, median absolute deviation of %d frames to estimate background model'%nframes,
                                  nframes*noffsets,
                                  parent,
                                  wx.PD_APP_MODAL|wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT|wx.PD_REMAINING_TIME)

        # allocate memory for median and mad
        self.med = num.zeros(framesize)
        self.mad = num.zeros(framesize)

        offseti = 0
        for imageoffset in range(0,nbytes,nbytessmall):
            if DEBUG: print "image offset = %d."%imageoffset

            # if this is the last iteration, there may not be npixelssmall left
            # store in npixelscurr the number of pixels that are to be read in,
            # store in seekperframecurr the amount to seek after reading.
            imageoffsetnext = imageoffset+nbytessmall
            if imageoffsetnext > nbytes:
                nbytescurr = nbytessmalllast
                seekperframecurr = seekperframelast
                imageoffsetnext = nbytes
            else:
                nbytescurr = nbytessmall
                seekperframecurr = seekperframe

            # allocate memory for buffer that holds data read in in current pass
            buf = num.zeros((nframes,nbytescurr),dtype=num.uint8)

            # seek to pixel imageoffset of the first frame
            fp.seek(imageoffset+headersize,0)

            if DEBUG: print 'Reading ...'

            # loop through frames
            for i in range(nframes):

                if params.interactive and not params.batch_executing:
                    (keepgoing,skip) = progressbar.Update(value=offseti*nframes + i,newmsg='Reading in piece of frame %d / %d, offset = %d (%d / %d)'%(i+1,nframes,imageoffset,offseti+1,noffsets))
                    if not keepgoing:
                        progressbar.Destroy()

                        if isbackup:
                            self.med = med0
                            self.mad = mad0
                        else:
                            delattr(self,'med')
                            delattr(self,'mad')

                        return False

                # seek to the desired part of the movie; 
                # skip bytesperchunk - the amount we read in in the last frame 
                if i > 0:
                    fp.seek(seekperframecurr,1)
                
                # read nrsmall rows
                data = fp.read(nbytescurr)
                if data == '':
                    raise NoMoreFramesException('EOF')
                buf[i,:] = num.fromstring(data,num.uint8)
                

            # compute the median and median absolute difference at each 
            # pixel location
            if DEBUG: print 'Computing ...'
            # sort all the histories to get the median
            buf = buf.transpose()
            buf.sort(axis=1,kind='mergesort')
            # store the median
            self.med[imageoffset:imageoffsetnext] = buf[:,middle1]
            if iseven:
                self.med[imageoffset:imageoffsetnext] += buf[:,middle2]
                self.med[imageoffset:imageoffsetnext] /= 2.
                
            # store the absolute difference
            buf = num.double(buf)
            for j in range(nframes):
                buf[:,j] = num.abs(buf[:,j] - self.med[imageoffset:imageoffsetnext])

            # sort
            buf.sort(axis=1,kind='mergesort')

            # store the median absolute difference
            self.mad[imageoffset:imageoffsetnext] = buf[:,middle1]
            if iseven:
                self.mad[imageoffset:imageoffsetnext] += buf[:,middle2]
                self.mad[imageoffset:imageoffsetnext] /= 2.

            #self.mad[imageoffset:imageoffsetnext] = buf[:,madorder1]*madweight1 + \
            #                                        buf[:,madorder2]*madweight2
            offseti += 1

        # estimate standard deviation assuming a Gaussian distribution
        # from the fact that half the data falls within mad
        # MADTOSTDFACTOR = 1./norminv(.75)
        MADTOSTDFACTOR = 1.482602
        self.mad *= MADTOSTDFACTOR
        
        self.mad.shape = params.movie_size
        self.med.shape = params.movie_size

        if params.interactive and not params.batch_executing:
            progressbar.Destroy()

        return True

    def est_bg( self, parent=None ):
        
        # make sure number of background frames is at most number of frames in movie
        if params.n_bg_frames > params.n_frames:
            params.n_bg_frames = params.n_frames

        # are we using the median or the mean?
        if params.use_median:
            if DEBUG: print 'params.movie.type = ' + str(params.movie.type)

            if params.movie.type == 'fmf' or \
                    (params.movie.type == 'avi' and \
                         params.movie.h_mov.bits_per_pixel == 8):
                succeeded = self.medmad(parent)
                if not succeeded:
                    return
                #(self.med,self.mad) = medmad(params.movie.filename,params.n_bg_frames,
                #                             params.movie_size[0],params.movie_size[1],
                #                             nframesskip,params.movie.h_mov.bytes_per_chunk,
                #                             params.movie.h_mov.chunk_start+8)
                self.center = self.med.copy()
                self.dev = self.mad.copy()
            elif params.movie.type == 'sbfmf':
                # the background model is fixed for this type of movie
                self.center = params.movie.h_mov.bgcenter.copy()
                self.dev = params.movie.h_mov.bgstd.copy()
                self.center.shape = params.movie.h_mov.framesize
                self.dev.shape = params.movie.h_mov.framesize
                if DEBUG: print 'center.shape = ' + str(self.center.shape)
                if DEBUG: print 'dev.shape = ' + str(self.dev.shape)
            else:
                succeeded = self.flexmedmad(parent)
                if not succeeded:
                    return
                self.center = self.med.copy()
                self.dev = self.mad.copy()
                #wx.MessageBox( "Only FMFs are supported for background estimation so far.",
                #               "Error", wx.ICON_ERROR )
        else:
            if params.movie.type == 'fmf' or params.movie.type == 'avi' \
                    or params.movie.type == 'avbin':
                succeeded = self.meanstd(parent)
                if not succeeded:
                    return
                #(self.mean,self.std) = meanstd(params.movie.filename,params.n_bg_frames,
                #                               params.movie_size[0],params.movie_size[1],
                #                               nframesskip,params.movie.h_mov.bytes_per_chunk,
                #                               params.movie.h_mov.chunk_start+4)
                self.center = self.mean.copy()
                self.dev = self.std.copy()
            elif params.movie.type == 'sbfmf':
                # the background model is fixed for this type of movie
                self.center = params.movie.h_mov.bgcenter.copy()
                self.dev = params.movie.h_mov.bgstd.copy()
                self.center.shape = params.movie.h_mov.framesize
                self.dev.shape = params.movie.h_mov.framesize
                if DEBUG: print 'center.shape = ' + str(self.center.shape)
                if DEBUG: print 'dev.shape = ' + str(self.dev.shape)
            else:
                wx.MessageBox( "Invalid movie type",
                               "Error", wx.ICON_ERROR )
        self.dev[self.dev < params.bg_std_min] = params.bg_std_min
        self.dev[self.dev > params.bg_std_max] = params.bg_std_max
        self.thresh = self.dev*params.n_bg_std_thresh
        self.thresh_low = self.dev*params.n_bg_std_thresh_low

        # image of which parts flies are not able to walk in
        self.UpdateIsArena()

        # approximate homomorphic filtering
        tmp = self.center.copy()
        issmall = tmp<1.
        tmp[issmall] = 1.
        self.hfnorm = self.hf.apply(self.center) / tmp
        self.hfnorm[issmall & (self.hfnorm<1.)] = 1.        

    def sub_bg(self,framenumber):
        """Reads image, subtracts background, then thresholds."""

        # read in the image
        # IndexError if requested frame < 0
        # NoMoreFramesException if frame is >= n_frames
        # ValueError if the frame size read in does not equal the buffer
        #   for uncompressed AVI
        
        self.curr_im, stamp = params.movie.get_frame( int(framenumber) )
        im = self.curr_im.astype( num.float )
        self.curr_stamp = stamp

        dfore = num.zeros(im.shape)
        if params.bg_type == params.BG_TYPE_LIGHTONDARK:
            # if white flies on a dark background, then we only care about differences
            # im - bgmodel.center > bgmodel.thresh
            dfore[self.isarena] = num.maximum(0,im[self.isarena] - self.center[self.isarena])
            
        elif params.bg_type == params.BG_TYPE_DARKONLIGHT:
            # if dark flies on a white background, then we only care about differences
            # bgmodel.center - im > bgmodel.thresh
            dfore[self.isarena] = num.maximum(0,self.center[self.isarena] - im[self.isarena])
            
        else:
            # otherwise, we care about both directions
            dfore[self.isarena] = num.abs(im[self.isarena] - self.center[self.isarena])

        dfore[self.isarena] /= self.dev[self.isarena]
        if params.n_bg_std_thresh > params.n_bg_std_thresh_low:
            bwhigh = num.zeros(im.shape,dtype=bool)
            bwlow = num.zeros(im.shape,dtype=bool)
            bwhigh[self.isarena] = dfore[self.isarena] > params.n_bg_std_thresh
            bwlow[self.isarena] = dfore[self.isarena] > params.n_bg_std_thresh_low
            
            # do hysteresis
            bw = morph.binary_propagation(bwhigh,mask=bwlow)
        else:
            bw = num.zeros(im.shape,dtype=bool)
            bw[self.isarena] = dfore[self.isarena] > params.n_bg_std_thresh

        # do morphology
        if params.do_use_morphology:
            if params.opening_radius > 0:
                bw = morph.binary_opening(bw,self.opening_struct)
            if params.closing_radius > 0:
                bw = morph.binary_closing(bw,self.closing_struct)

        return (dfore,bw)

    def create_morph_struct(self, radius):
        if radius <= 0:
            s = False
        elif radius == 1:
            s = num.ones((2,2),bool)
        else:
            s = morph.generate_binary_structure(2,1)
            s = morph.iterate_structure(s,radius)
        return s

    def ShowBG( self, parent, framenumber, old_thresh ):
        
        self.old_thresh = params.n_bg_std_thresh
        self.show_frame = framenumber

        rsrc = xrc.XmlResource( THRESH_RSRC_FILE )
        self.frame = rsrc.LoadFrame( parent, "frame_Ctrax_bg" )

        # calculate background image, if not already done
        if not hasattr( self, 'center' ):
            wx.BeginBusyCursor()
            wx.Yield()
            self.est_bg(self.frame)
            wx.EndBusyCursor()

        self.hf.SetFrame(self.frame)
        self.hf.frame.Bind( wx.EVT_BUTTON, self.OnQuitHF, id=xrc.XRCID("hm_done_button") )
        self.hf.frame.Bind( wx.EVT_CLOSE, self.OnQuitHF )

        # event bindings
        self.menu = self.frame.GetMenuBar()
        #self.frame.Bind( wx.EVT_MENU, self.OnShowBG, id=xrc.XRCID("menu_show_bg") )
        #self.frame.Bind( wx.EVT_MENU, self.OnShowThresh, id=xrc.XRCID("menu_show_thresh") )
        #self.frame.Bind( wx.EVT_MENU, self.OnFrameSlider, id=xrc.XRCID("menu_show_bin") )
        #self.frame.Bind( wx.EVT_BUTTON, self.OnCalcButton, id=xrc.XRCID("button_calculate") )
        #self.frame.Bind( wx.EVT_BUTTON, self.OnResetButton, id=xrc.XRCID("button_reset") )
        
        # control handles
        self.frame_text = xrc.XRCCTRL( self.frame, "text_frame" )
        #self.reset_button = xrc.XRCCTRL( self.frame, "button_reset" )
        #if self.old_thresh is None:
        #    self.reset_button.Enable( False )
        
        # set threshold
        self.thresh_slider = xrc.XRCCTRL( self.frame, "slider_thresh" )
        self.thresh_slider.SetScrollbar(self.GetThresholdScrollbar(),0,params.sliderres,25)
        self.thresh_low_slider = xrc.XRCCTRL( self.frame, "slider_low_thresh" )
        self.thresh_low_slider.SetScrollbar(self.GetThresholdLowScrollbar(),0,params.sliderres,25)

        # set the bounds of the slider
        self.SetThreshBounds()

        # we can also input the threshold by typing it; make sure slider and text input
        # show the same value
        self.thresh_textinput = xrc.XRCCTRL( self.frame, "threshold_text_input" )
        self.frame.Bind( wx.EVT_SCROLL, self.OnThreshSlider, self.thresh_slider )
        wxvt.setup_validated_float_callback( self.thresh_textinput,
                                             xrc.XRCID("threshold_text_input"),
                                             self.OnThreshTextEnter,
                                             pending_color=params.wxvt_bg )
        self.thresh_low_textinput = xrc.XRCCTRL( self.frame, "low_threshold_text_input" )
        self.frame.Bind( wx.EVT_SCROLL, self.OnThreshSlider, self.thresh_low_slider )
        wxvt.setup_validated_float_callback( self.thresh_low_textinput,
                                             xrc.XRCID("low_threshold_text_input"),
                                             self.OnThreshTextEnter,
                                             pending_color=params.wxvt_bg )
        self.SetThreshold()
        self.frame_slider = xrc.XRCCTRL( self.frame, "slider_frame" )
        self.frame.Bind( wx.EVT_SCROLL, self.OnFrameSlider, self.frame_slider )
        self.frame_slider.SetScrollbar(self.show_frame,0,params.n_frames-1,params.n_frames/10)
        # which image are we showing
        self.bg_img_chooser = xrc.XRCCTRL( self.frame, "view_type_input" )
        self.bg_img_chooser.SetSelection(self.show_img_type)
        self.frame.Bind( wx.EVT_CHOICE, self.OnImgChoice, self.bg_img_chooser)
        # background type
        self.bg_type_chooser = xrc.XRCCTRL(self.frame,"bg_type_input")
        self.bg_type_chooser.SetSelection(params.bg_type)
        self.frame.Bind(wx.EVT_CHOICE,self.OnBgTypeChoice, self.bg_type_chooser)
        # minimum bg std
        self.bg_minstd_textinput = xrc.XRCCTRL(self.frame,"bg_min_std_input")
        wxvt.setup_validated_float_callback( self.bg_minstd_textinput,
                                             xrc.XRCID("bg_min_std_input"),
                                             self.OnMinStdTextEnter,
                                             pending_color=params.wxvt_bg )
        self.bg_minstd_textinput.SetValue( '%.2f'%params.bg_std_min )
        # maximum bg std
        self.bg_maxstd_textinput = xrc.XRCCTRL(self.frame,"bg_max_std_input")
        wxvt.setup_validated_float_callback( self.bg_maxstd_textinput,
                                             xrc.XRCID("bg_max_std_input"),
                                             self.OnMinStdTextEnter,
                                             pending_color=params.wxvt_bg )
        self.bg_maxstd_textinput.SetValue( '%.2f'%params.bg_std_max )
        # normalization type
        self.bg_norm_chooser = xrc.XRCCTRL(self.frame,"bg_normalization_input")
        self.bg_norm_chooser.SetSelection(params.bg_norm_type)
        self.frame.Bind( wx.EVT_CHOICE, self.OnNormChoice, self.bg_norm_chooser)
        # homomorphic filter settings dialog
        self.hf_button = xrc.XRCCTRL(self.frame,"homomorphic_settings")
        self.hf_button.Enable(False)
        self.frame.Bind(wx.EVT_BUTTON,self.OnHFClick,self.hf_button)
        self.hf_window_open = False

        # detect arena button
        self.detect_arena_button = xrc.XRCCTRL(self.frame,"detect_arena_button")
        self.detect_arena_checkbox = xrc.XRCCTRL(self.frame,"detect_arena_checkbox")
        self.detect_arena_button.Enable(params.do_set_circular_arena)
        self.detect_arena_checkbox.SetValue(params.do_set_circular_arena)
        
        self.frame.Bind(wx.EVT_BUTTON,self.OnDetectArenaClick,self.detect_arena_button)
        self.frame.Bind(wx.EVT_CHECKBOX,self.OnDetectArenaCheck,self.detect_arena_checkbox)
        self.detect_arena_window_open = False
        
        # min value for isarena
        self.min_nonarena_textinput = xrc.XRCCTRL(self.frame,"min_nonfore_intensity_input")
        wxvt.setup_validated_float_callback( self.min_nonarena_textinput,
                                             xrc.XRCID("min_nonfore_intensity_input"),
                                             self.OnMinNonArenaTextEnter,
                                             pending_color=params.wxvt_bg )
        self.min_nonarena_textinput.SetValue( '%.1f'%params.min_nonarena )
        # min value for isarena
        self.max_nonarena_textinput = xrc.XRCCTRL(self.frame,"max_nonfore_intensity_input")
        wxvt.setup_validated_float_callback( self.max_nonarena_textinput,
                                             xrc.XRCID("max_nonfore_intensity_input"),
                                             self.OnMaxNonArenaTextEnter,
                                             pending_color=params.wxvt_bg )
        self.max_nonarena_textinput.SetValue( '%.1f'%params.max_nonarena )

        # morphology
        self.morphology_checkbox = xrc.XRCCTRL(self.frame,"morphology_checkbox")
        self.morphology_checkbox.SetValue(params.do_use_morphology)
        self.frame.Bind(wx.EVT_CHECKBOX,self.OnMorphologyCheck,self.morphology_checkbox)
        self.opening_radius_textinput = xrc.XRCCTRL(self.frame,"opening_radius")
        wxvt.setup_validated_integer_callback( self.opening_radius_textinput,
                                               xrc.XRCID("opening_radius"),
                                               self.OnOpeningRadiusTextEnter,
                                               pending_color=params.wxvt_bg )
        self.opening_radius_textinput.SetValue( '%d'%params.opening_radius )
        self.opening_radius_textinput.Enable(params.do_use_morphology)
        self.closing_radius_textinput = xrc.XRCCTRL(self.frame,"closing_radius")
        self.opening_struct = self.create_morph_struct(params.opening_radius)
        wxvt.setup_validated_integer_callback( self.closing_radius_textinput,
                                               xrc.XRCID("closing_radius"),
                                             self.OnClosingRadiusTextEnter,
                                             pending_color=params.wxvt_bg )
        self.closing_radius_textinput.SetValue( '%d'%params.closing_radius )
        self.closing_radius_textinput.Enable(params.do_use_morphology)
        self.closing_struct = self.create_morph_struct(params.closing_radius)

        self.fixbg_button = xrc.XRCCTRL(self.frame,"fixbg_button")
        self.frame.Bind(wx.EVT_BUTTON,self.OnFixBgClick,self.fixbg_button)
        self.fixbg_window_open = False

        self.roi_button = xrc.XRCCTRL(self.frame,"roi_button")
        self.frame.Bind(wx.EVT_BUTTON,self.OnROIClick,self.roi_button)
        self.roi_window_open = False

        # make image window
        self.img_panel = xrc.XRCCTRL( self.frame, "panel_img" )
        box = wx.BoxSizer( wx.VERTICAL )
        self.img_panel.SetSizer( box )
        self.img_wind = wxvideo.DynamicImageCanvas( self.img_panel, -1 )
        self.img_wind.set_resize(True)
        box.Add( self.img_wind, 1, wx.EXPAND )
        self.img_panel.SetAutoLayout( True )
        self.img_panel.Layout()
        
    def SetThreshBounds( self , maxv=None):
        if maxv is not None:
            # changing the conversion from scrollbar units to threshold units
            self.scrollbar2thresh = float(maxv) / float(params.sliderres)
            return
        (self.dfore,self.bw) = self.sub_bg(self.show_frame)
        # what is the maximum dfore?
        self.max_slider = num.max(self.dfore)
        self.scrollbar2thresh = float(self.max_slider) / float(params.sliderres)

    def SetThreshold(self):
        #self.thresh_textinput.SetValue( '%.1f'%params.n_bg_std_thresh )
        if params.n_bg_std_thresh < params.n_bg_std_thresh_low:
            params.n_bg_std_thresh_low = params.n_bg_std_thresh
        self.thresh_slider.SetThumbPosition( self.GetThresholdScrollbar() )
        self.thresh_low_slider.SetThumbPosition( self.GetThresholdLowScrollbar() )
        self.thresh_textinput.SetValue('%.1f'%params.n_bg_std_thresh)
        self.thresh_low_textinput.SetValue('%.1f'%params.n_bg_std_thresh_low)
    def GetThresholdScrollbar(self):
        if not hasattr(self,'scrollbar2thresh'):
            self.SetThreshBounds()
        if self.scrollbar2thresh == 0:
            return 0
        else:
            return(params.n_bg_std_thresh/self.scrollbar2thresh)

    def GetThresholdLowScrollbar(self):
        if not hasattr(self,'scrollbar2thresh'):
            self.SetThreshBounds()
        if self.scrollbar2thresh == 0:
            return 0
        else:
            return(params.n_bg_std_thresh_low/self.scrollbar2thresh)

    def ReadScrollbar(self):
        return (self.thresh_slider.GetThumbPosition()*self.scrollbar2thresh)


    def ReadScrollbarLow(self):
        return (self.thresh_low_slider.GetThumbPosition()*self.scrollbar2thresh)

    #def OnCalcButton( self, evt ):
    #    wx.BeginBusyCursor()
    #    wx.Yield()
    #    #self.est_thresh()
    #    wx.EndBusyCursor()
    #    wx.Yield()
    #    #
    #    #self.thresh_slider.SetValue( int(round(self.bg_thresh*255.)) )
    #    #self.OnThreshSlider( True )

    #def OnResetButton( self, evt ):
    #    self.thresh_slider.SetValue( self.old_thresh*params.sliderres )
    #    self.OnThreshSlider( True )

    def OnShowBG( self, evt ):
        #self.menu.Enable( xrc.XRCID("menu_show_bin"), False )
        self.frame_slider.Enable( False )
        self.DoSub()
        
    def OnShowThresh( self, evt ):
        #self.menu.Enable( xrc.XRCID("menu_show_bin"), True )
        self.frame_slider.Enable( True )
        self.DoSub()

    def OnImgChoice( self, evt ):
        new_img_type = self.bg_img_chooser.GetSelection()
        if new_img_type is not wx.NOT_FOUND:
            self.show_img_type = new_img_type
            self.DoSub()

    def OnBgTypeChoice( self, evt ):
        new_bg_type = self.bg_type_chooser.GetSelection()
        if new_bg_type is not wx.NOT_FOUND:
            params.bg_type = new_bg_type
            self.DoSub()
            # change threshold bounds
            self.SetThreshBounds()
            if params.n_bg_std_thresh > self.max_slider:
                params.n_bg_std_thresh = self.max_slider
                self.SetThreshold()
                self.DoSub()

    def OnNormChoice( self, evt ):
        new_norm_type = self.bg_norm_chooser.GetSelection()
        if new_norm_type is not wx.NOT_FOUND:
            params.bg_norm_type = new_norm_type
            if evt is not None:
                if params.bg_norm_type == params.BG_NORM_STD:
                    if params.use_median:
                        self.dev[:] = self.mad
                        self.dev[self.dev < params.bg_std_min] = params.bg_std_min
                        self.dev[self.dev > params.bg_std_max] = params.bg_std_max
                    else:
                        self.dev[:] = self.std
                        self.dev[self.dev < params.bg_std_min] = params.bg_std_min
                        self.dev[self.dev > params.bg_std_max] = params.bg_std_max
                    self.hf_button.Enable(False)
                elif params.bg_norm_type == params.BG_NORM_INTENSITY:
                    self.dev[:] = self.center
                    self.hf_button.Enable(False)
                elif params.bg_norm_type == params.BG_NORM_HOMOMORPHIC:
                    self.dev[:] = self.hfnorm
                    self.hf_button.Enable(True)
                # get the old maxdfore
                if hasattr(self,'dfore'):
                    oldmaxdfore = num.max(self.dfore)
                    donormalize = True
                else:
                    donormalize = False
                (self.dfore,self.bw) = self.sub_bg(self.show_frame)
                # change the scale of the threshold
                self.SetThreshBounds()
                if donormalize:
                    params.n_bg_std_thresh *= self.max_slider/oldmaxdfore
                    params.n_bg_std_thresh_low *= self.max_slider/oldmaxdfore
                self.SetThreshold()
                self.DoSub()

    def OnHFClick( self, evt ):
        if self.hf_window_open:
            self.hf.frame.Raise()
            return
        self.hf_window_open = True
        self.hf.frame.Show()

    def OnDetectArenaCheck( self, evt ):
        if evt is None:
            return
        params.do_set_circular_arena = self.detect_arena_checkbox.GetValue()
        self.detect_arena_button.Enable(params.do_set_circular_arena)
        self.UpdateIsArena()
        self.DoSub()

    def OnMorphologyCheck( self, evt ):
        if evt is None:
            return
        params.do_use_morphology = self.morphology_checkbox.GetValue()
        self.opening_radius_textinput.Enable(params.do_use_morphology)
        self.closing_radius_textinput.Enable(params.do_use_morphology)
        self.DoSub()

    def OnOpeningRadiusTextEnter(self, evt):
        textentered = self.opening_radius_textinput.GetValue()
        # convert to number
        try:
            params.opening_radius = int(textentered)
        except ValueError:
            # if not a number, then set the value to the previous threshold
            self.opening_radius_textinput.SetValue('%d'%params.opening_radius)
            return
        self.opening_struct = self.create_morph_struct(params.opening_radius)
        if evt is not None:
            self.DoSub()

    def OnClosingRadiusTextEnter(self, evt):
        textentered = self.closing_radius_textinput.GetValue()
        # convert to number
        try:
            params.closing_radius = int(textentered)
        except ValueError:
            # if not a number, then set the value to the previous threshold
            self.closing_radius_textinput.SetValue('%d'%params.closing_radius)
            return
        self.closing_struct = self.create_morph_struct(params.closing_radius)

        if evt is not None:
            self.DoSub()

    def UpdateIsArena(self):
        # update isarena
        self.isarena = num.ones(params.movie_size,num.bool)
        if hasattr(params,'max_nonarena'):
            self.isarena = self.center > params.max_nonarena
        if hasattr(params,'min_nonarena'):
            self.isarena = self.isarena & \
                              (self.center < params.min_nonarena)
        if params.do_set_circular_arena and \
           hasattr(params,'arena_center_x') and \
           (params.arena_center_x is not None):
            self.isarena = self.isarena & ( ((params.GRID.X - params.arena_center_x)**2. + (params.GRID.Y - params.arena_center_y)**2) <= params.arena_radius**2. )
        if hasattr(self,'roi'):
            self.isarena = self.isarena & self.roi

    def OnDetectArenaClick( self, evt ):
        if self.detect_arena_window_open:
            self.detect_arena.frame.Raise()
            return
        self.detect_arena_window_open = True
        self.detect_arena = setarena.SetArena(self.frame,self.center)
        self.detect_arena.frame.Bind( wx.EVT_BUTTON, self.OnQuitDetectArena, id=xrc.XRCID("done_button") )
        self.detect_arena.frame.Bind( wx.EVT_CLOSE, self.OnQuitDetectArena )
        self.detect_arena.frame.Show()

    def OnQuitDetectArena( self, evt ):

        # save parameters found
        [params.arena_center_x,params.arena_center_y,
         params.arena_radius] = self.detect_arena.GetArenaParameters()

        # update isarena
        self.UpdateIsArena()

        # close frame
        self.detect_arena_window_open = False
        self.detect_arena.frame.Destroy()
        delattr(self,'detect_arena')

        # redraw
        self.DoSub()

    def OnFixBgClick( self, evt ):
        if self.fixbg_window_open:
            self.fixbg.frame.Raise()
            return
        self.fixbg_window_open = True
        self.fixbg = fixbg.FixBG(self.frame,self)
        self.fixbg.frame.Bind( wx.EVT_BUTTON, self.OnQuitFixBg, id=xrc.XRCID("quit_button") )
        self.fixbg.frame.Bind( wx.EVT_CLOSE, self.OnQuitFixBg )
        self.fixbg.frame.Show()

    def OnQuitFixBg( self, evt ):

        # close frame
        reallyquit = self.fixbg.CheckSave()

        if not reallyquit:
            return

        self.fixbg_window_open = False
        self.fixbg.frame.Destroy()
        delattr(self,'fixbg')

        # redraw
        self.DoSub()

    def SetCenter(self,newcenter):
        self.center = newcenter.copy()
        if params.use_median:
            self.med = newcenter.copy()
        else:
            self.mean = newcenter.copy()

    def SetDev(self,newdev):
        self.dev = newdev.copy()
        if params.use_median:
            self.mad = newdev.copy()
        else:
            self.std = newdev.copy()

    def SetFixBgPolygons(self,polygons):
        self.fixbg_polygons = polygons[:]

    def OnROIClick( self, evt ):
        if self.roi_window_open:
            self.roigui.frame.Raise()
            return
        self.roi_window_open = True
        self.roigui = roi.ROI(self.frame,self)
        self.roigui.frame.Bind( wx.EVT_BUTTON, self.OnQuitROI, id=xrc.XRCID("quit_button") )
        self.roigui.frame.Bind( wx.EVT_CLOSE, self.OnQuitROI )
        self.roigui.frame.Show()

    def OnQuitROI( self, evt ):

        # close frame
        reallyquit = self.roigui.CheckSave()

        if not reallyquit:
            return

        self.roi_window_open = False
        self.roigui.frame.Destroy()
        delattr(self,'roigui')

        # redraw
        self.DoSub()

    def SetROI(self,roi):
        self.roi = roi.copy()
        self.UpdateIsArena()

    def OnQuitHF(self, evt ):
        self.hf_window_open = False
        self.hf.frame.Hide()
        self.DoSub()
    
    def OnThreshSlider( self, evt ):
        params.n_bg_std_thresh = self.ReadScrollbar()
        params.n_bg_std_thresh_low = self.ReadScrollbarLow()
        if params.n_bg_std_thresh_low > params.n_bg_std_thresh:
            params.n_bg_std_thresh_low = params.n_bg_std_thresh
        self.thresh_textinput.SetValue( '%.1f'%params.n_bg_std_thresh )
        self.thresh_low_textinput.SetValue( '%.1f'%params.n_bg_std_thresh_low )
        if evt is not None: self.DoSub()

    def OnThreshTextEnter( self, evt ):
        # get the value entered
        params.n_bg_std_thresh = float(self.thresh_textinput.GetValue())
        params.n_bg_std_thresh_low = float(self.thresh_low_textinput.GetValue())

        # make sure it is in a valid interval
        if params.n_bg_std_thresh < 0: params.n_bg_std_thresh = 0
        if params.n_bg_std_thresh_low < 0: params.n_bg_std_thresh_low = 0
        if params.n_bg_std_thresh_low > params.n_bg_std_thresh:
            params.n_bg_std_thresh_low = params.n_bg_std_thresh

        self.thresh_slider.SetThumbPosition( self.GetThresholdScrollbar() )
        self.thresh_low_slider.SetThumbPosition( self.GetThresholdLowScrollbar() )

        if evt is not None: self.DoSub()

    def OnMinStdTextEnter( self, evt ):
        # get the value entered
        textentered = self.bg_minstd_textinput.GetValue()
        maxtextentered = self.bg_maxstd_textinput.GetValue()
        # convert to number
        try:
            numentered = float(textentered)
            maxnumentered = float(maxtextentered)

            # make sure that min <= max
            if numentered > maxnumentered:
                self.bg_minstd_textinput.SetValue('%.2f'%params.bg_std_min)
                self.bg_maxstd_textinput.SetValue('%.2f'%params.bg_std_max)
                return

            # make sure it is in a valid interval
            if numentered < 0:
                numentered = 0
            if maxnumentered < 0:
                maxnumentered = 0
            if numentered > 255:
                maxnumentered = 255
            params.bg_std_min = numentered
            params.bg_std_max = maxnumentered
            self.bg_minstd_textinput.SetValue('%.2f'%params.bg_std_min)
            self.bg_maxstd_textinput.SetValue('%.2f'%params.bg_std_max)
            if params.use_median:
                self.dev[:] = self.mad
            else:
                self.dev[:] = self.std
            self.dev[self.dev < params.bg_std_min] = params.bg_std_min
            self.dev[self.dev > params.bg_std_max] = params.bg_std_max
        except ValueError:
            # if not a number, then set the value to the previous threshold
            self.bg_minstd_textinput.SetValue('%.2f'%params.bg_std_min)
            self.bg_maxstd_textinput.SetValue('%.2f'%params.bg_std_max)
            return
        if evt is not None:
            self.DoSub()
            self.SetThreshBounds()
            if params.n_bg_std_thresh > self.max_slider:
                params.n_bg_std_thresh = self.max_slider
                self.SetThreshold()
                self.DoSub()

    def OnMinNonArenaTextEnter( self, evt ):
        # get the value entered
        textentered = self.min_nonarena_textinput.GetValue()
        # convert to number
        try:
            params.min_nonarena = float(textentered)
        except ValueError:
            # if not a number, then set the value to the previous threshold
            self.min_nonarena_textinput.SetValue('%.1f'%params.min_nonarena)
            return
        if evt is not None:
            self.UpdateIsArena()
            self.DoSub()

    def OnMaxNonArenaTextEnter( self, evt ):
        # get the value entered
        textentered = self.max_nonarena_textinput.GetValue()
        # convert to number
        try:
            params.max_nonarena = float(textentered)
        except ValueError:
            # if not a number, then set the value to the previous threshold
            self.max_nonarena_textinput.SetValue('%.1f'%params.max_nonarena)
            return
        if evt is not None:
            self.UpdateIsArena()
            self.DoSub()


    def OnFrameSlider( self, evt ):
        self.show_frame = self.frame_slider.GetThumbPosition()
        self.frame_text.SetLabel( "frame %d"%(self.show_frame) )
        if evt is not None:
            self.DoSub()
            self.SetThreshBounds()

    def DoSub( self ):
        """Do background subtraction."""
        wx.BeginBusyCursor()
        wx.Yield()
        (self.dfore,self.bw) = self.sub_bg(self.show_frame)
        self.windowsize = [self.img_panel.GetRect().GetHeight(),self.img_panel.GetRect().GetWidth()]
        if self.show_img_type == params.SHOW_BACKGROUND:
            img_8 = imagesk.double2mono8(self.center,donormalize=False)
            self.img_wind.update_image_and_drawings('bg',img_8,format='MONO8')
        elif self.show_img_type == params.SHOW_DISTANCE:
            #self.img = self.dfore.copy()
            img_8 = imagesk.double2mono8(self.dfore)
            self.img_wind.update_image_and_drawings('bg',img_8,format='MONO8')
        elif self.show_img_type == params.SHOW_THRESH:
            img_8 = imagesk.double2mono8(self.bw.astype(float))
            self.img_wind.update_image_and_drawings('bg',img_8,format='MONO8')
        elif self.show_img_type == params.SHOW_NONFORE:
            isnotarena = self.isarena == False
            img_8 = imagesk.double2mono8(isnotarena.astype(float))
            self.img_wind.update_image_and_drawings('bg',img_8,format='MONO8')
        elif self.show_img_type == params.SHOW_DEV:
            mu = num.mean(self.dev)
            sig = num.std(self.dev)
            if sig == 0:
                resp = wx.MessageBox("Normalization image is uniformly " + str(mu),"Uniform Normalization Image",wx.OK)
                img_8 = imagesk.double2mono8(self.dev)
            else:
                n1 = mu - 3.*sig
                n2 = mu + 3.*sig
                img_8 = imagesk.double2mono8(self.dev.clip(n1,n2))
            self.img_wind.update_image_and_drawings('bg',img_8,format='MONO8')
        elif self.show_img_type == params.SHOW_CC:
            (L,ncc) = meas.label(self.bw)
            #self.img = colormap_image(L)
            img_8 = colormap_image(L)
            self.img_wind.update_image_and_drawings('bg',img_8,format='RGB8')
        elif self.show_img_type == params.SHOW_ELLIPSES:
            # find ellipse observations
            obs = find_ellipses(self.dfore,self.bw,False)
            im, stamp = params.movie.get_frame( int(self.show_frame) )
            img_8 = imagesk.double2mono8(im,donormalize=False)
            plot_linesegs = draw_ellipses(obs)
            (linesegs,linecolors) = \
                imagesk.separate_linesegs_colors(plot_linesegs)
            img_8 = imagesk.double2mono8(im,donormalize=False)
            self.img_wind.update_image_and_drawings('bg',
                                                    img_8,
                                                    format='MONO8',
                                                    linesegs=linesegs,
                                                    lineseg_colors=linecolors
                                                    )

        self.img_wind.Refresh(eraseBackground=False)
        wx.EndBusyCursor()

class BgSettingsDialog:

    def __init__(self,parent,bg):
        #self.const = bg.const
        #const.n_frames = bg.n_frames
        self.bg = bg
        rsrc = xrc.XmlResource( SETTINGS_RSRC_FILE )
        self.frame = rsrc.LoadFrame(parent, "bg_settings")
        self.algorithm = xrc.XRCCTRL( self.frame, "algorithm" )
        self.nframes = xrc.XRCCTRL( self.frame, "nframes" )
        self.bg_firstframe = xrc.XRCCTRL( self.frame, "bg_firstframe" )
        self.bg_lastframe = xrc.XRCCTRL( self.frame, "bg_lastframe" )
        self.calculate = xrc.XRCCTRL( self.frame, "calculate_button" )
        if params.use_median:
            self.algorithm.SetSelection(params.ALGORITHM_MEDIAN)
        else:
            self.algorithm.SetSelection(params.ALGORITHM_MEAN)
        self.nframes.SetValue( str(params.n_bg_frames))
        wxvt.setup_validated_integer_callback( self.nframes,
                                               xrc.XRCID("nframes"),
                                               self.OnTextValidated,
                                               pending_color=params.wxvt_bg )
        self.bg_firstframe.SetValue( str(params.bg_firstframe))
        wxvt.setup_validated_integer_callback( self.bg_firstframe,
                                               xrc.XRCID("bg_firstframe"),
                                               self.OnTextValidatedFirstframe,
                                               pending_color=params.wxvt_bg )
        self.bg_lastframe.SetValue( str(params.bg_lastframe))
        wxvt.setup_validated_integer_callback( self.bg_lastframe,
                                               xrc.XRCID("bg_lastframe"),
                                               self.OnTextValidatedLastframe,
                                               pending_color=params.wxvt_bg )
        self.frame.Bind( wx.EVT_CHOICE, self.OnAlgorithmChoice, self.algorithm)
        self.frame.Bind( wx.EVT_BUTTON, self.OnCalculate, self.calculate )
        self.frame.Show()

    def OnCalculate(self,evt):
        wx.BeginBusyCursor()
        wx.Yield()
        self.bg.est_bg(self.frame)
        wx.EndBusyCursor()

    def OnTextValidated(self,evt):
        tmp = int(self.nframes.GetValue())
        if tmp > 0:
            params.n_bg_frames = tmp
            if params.n_bg_frames > params.n_frames:
                params.n_bg_frames = params.n_frames
            self.nframes.SetValue(str(params.n_bg_frames))

    def OnTextValidatedFirstframe(self,evt):
        tmp0 = int(self.bg_firstframe.GetValue())
        tmp = tmp0
        tmp = max(tmp,0)
        tmp = min(tmp,params.bg_lastframe)
        tmp = min(tmp,params.n_frames-1)
        params.bg_firstframe = tmp
        if not tmp == tmp0:
            self.bg_firstframe.SetValue(str(params.bg_firstframe))

    def OnTextValidatedLastframe(self,evt):
        tmp0 = int(self.bg_lastframe.GetValue())
        tmp = tmp0
        tmp = max(tmp,params.bg_firstframe)
        params.bg_lastframe = tmp
        if not tmp == tmp0:
            self.bg_lastframe.SetValue(str(params.bg_lastframe))

    def OnAlgorithmChoice(self,evt):
        new_alg_type = self.algorithm.GetSelection()
        if new_alg_type is not wx.NOT_FOUND:
            if new_alg_type == params.ALGORITHM_MEDIAN:
                params.use_median = True
            else:
                params.use_median = False

#class StoreThread (threading.Thread):
#    """This thread's job is to read images, homomorphically filter them,
#    and place them in memory at a specified location."""
#    def __init__( self, hf, framelist, raw_img_mem, img_mem, indices, lock ):
#        threading.Thread.__init__( self )
#        self.filt = hf
#        self.framelist = framelist
#        self.img_mem = img_mem
#        self.raw_img_mem = raw_img_mem
#        self.indices = indices
#        self.lock = lock
#
#    def run( self ):
#        self.lock.acquire()
#        try:
#            for ff in range( len(self.framelist) ):
#                if self.raw_img_mem is not None:
#                    raw, filt = self.filt.filt_and_raw( frame )
#                    self.img_mem[self.indices[ff],:,:] = filt
#                    self.raw_img_mem[self.indices[ff],:,:] = raw
#                else:
#                    self.img_mem[:,:,self.indices[ff]] = \
#                                self.filt.homo_filt( self.framelist[ff] )
#        finally:
#            self.lock.release()
#
#class MaxMinThread (threading.Thread):
#    """This thread adds each filtered image to an accumulator and then
#    also sets global max and min values in specified arrays."""
#    def __init__( self, hf, framelist, lock,
#                  raw_acc, raw_acc_lock, filt_acc, filt_acc_lock,
#                  maxv, maxv_lock, minv, minv_lock ):
#        threading.Thread.__init__( self )
#        self.filt = hf
#        self.framelist = framelist
#        self.lock = lock
#
#        self.raw_acc = raw_acc
#        self.raw_acc_lock = raw_acc_lock
#        self.filt_acc = filt_acc
#        self.filt_acc_lock = filt_acc_lock
#        self.maxv = maxv
#        self.maxv_lock = maxv_lock
#        self.minv = minv
#        self.minv_lock = minv_lock
#
#    def run( self ):
#        self.lock.acquire()
#        try:
#            for frame in self.framelist:
#                # get image
#                raw, filt = self.filt.filt_and_raw_float( frame )
#                # add to accumulators
#                self.raw_acc_lock.acquire()
#                self.raw_acc += raw
#                self.raw_acc_lock.release()
#                self.filt_acc_lock.acquire()
#                self.filt_acc += filt
#                self.filt_acc_lock.release()
#                # compare to max-val array
#                f = filt > self.maxv
#                self.maxv_lock.acquire()
#                self.maxv[f] = filt[f]
#                self.maxv_lock.release()
#                # compare to min-val array
#                f = filt < self.minv
#                self.minv_lock.acquire()
#                self.minv[f] = filt[f]
#                self.minv_lock.release()
#        finally:
#            self.lock.release()
#
#class HistThread (threading.Thread):
#    """This thread creates histograms of the error between a background image
#    and a new image."""
#    def __init__( self, hf, framelist, lock, bg_img, binedges, n_bins,
#                  bincount, bincount_lock ):
#        threading.Thread.__init__( self )
#        self.filt = hf
#        self.framelist = framelist
#        self.lock = lock
#        
#        self.bg_img = bg_img
#        self.binedges = binedges
#        self.n_bins = n_bins
#
#        self.bincount = bincount
#        self.bincount_lock = bincount_lock
#
#    def run( self ):
#        self.lock.acquire()
#        N = self.filt.movie.get_width() * self.filt.movie.get_height()
#        try:
#            for frame in self.framelist:
#                # get image
#                img = self.filt.homo_filt( frame )
#                # compute error
#                err = num.abs( img - self.bg_img )
#                # choose NHISTPERFRAME randomly
#                r = num.ceil( N * random.rand( self.n_bins ) ).astype( num.int )
#                tmp, bins = num.histogram( err.flat[r], self.binedges )
#                # add to count
#                self.bincount_lock.acquire()
#                self.bincount += tmp[:-1]
#                self.bincount_lock.release()
#                # choose NHISTPERFRAME with highest error
#                sortederr = num.sort( make_vect( err ) )
#                desc_sortederr = sortederr[-1::-1]
#                tmp, bins = num.histogram( desc_sortederr[:self.n_bins], self.binedges )
#                self.bincount_lock.acquire()
#                self.bincount += tmp[:-1]
#                self.bincount_lock.release()
#        finally:
#            self.lock.release()

#def save_image( filename, img ):
#    """Saves an image as a MATLAB array."""
#    scipy.io.savemat( filename + '.mat', {filename: img})
#    print "saved", filename, ".mat"

#def read_image( filename ):
#    """Reads an image from a mat-file."""
#    new_dict = scipy.io.loadmat( filename )
#    print "loaded", filename, ".mat"
#    return new_dict[filename]

#def write_img_text( img, filename ):
#    fp = open( filename, "w" )
#    for r in img:
#        for c in r:
#            fp.write( "%f\t"%c )
#        fp.write( "\n" )
#    fp.close()

