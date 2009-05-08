# algorithm.py
# KMB 09/07/07

import numpy as num
from scipy.linalg.basic import eps
import time
import wx
from wx import xrc
import setarena
import chooseorientations
import os

import wxvalidatedtext as wxvt # part of Motmot

from version import DEBUG
import annfiles as annot
import ellipsesk as ell
import settings
import hindsight
import sys

from params import params

# import pkg_resources # part of setuptools
SETTINGS_RSRC_FILE = "track_settings.xrc"

# CtraxApp.Track #################################################
class CtraxAlgorithm (settings.AppWithSettings):
    """Cannot be used alone -- this class only exists
    to keep algorithm code together in one file."""
    def Track( self ):
        """Run the m-tracker."""

        ## initialization ##

        # maximum number of frames we will look back to fix errors
        self.maxlookback = max(params.lostdetection_length,
                               params.spuriousdetection_length,
                               params.mergeddetection_length,
                               params.splitdetection_length)

        # initialize annotation file
        self.ann_file = annot.AnnotationFile( self.ann_filename, params.version, True )
        self.ann_file.WriteAnnHeader( params.start_frame, self.bg_imgs )

        # if we haven't done any tracking yet
        if (self.ann_data is None) or (len(self.ann_data) == 0):
            self.ann_data = []
            params.nids = 0

        # write the current tracks to file
        for ff in range(len(self.ann_data)-self.maxlookback):
            self.ann_file.write_ellipses(self.ann_data[ff])
        self.lastframewritten = max(-1,len(self.ann_data)-self.maxlookback-1)

        wx.Yield()

        # initialize hindsight data structures
        self.hindsight = hindsight.Hindsight(self.ann_data,self.bg_imgs)
        self.hindsight.initialize_milestones()

        self.break_flag = False

        #print 'In Track, just before main loop: self.start_frame = %d, params.start_frame = %d, len(ann_data) = %d'%(self.start_frame,params.start_frame,len(self.ann_data))

        for self.start_frame in range(self.start_frame,self.n_frames):
        
            if self.break_flag:
                break

            last_time = time.time()

            # perform background subtraction
            try:
                (self.dfore,self.isfore) = self.bg_imgs.sub_bg( self.start_frame )
            except:
                break

            # write to sbfmf
            if self.dowritesbfmf:
                self.movie.writesbfmf_writeframe(self.isfore,
                                                 self.bg_imgs.curr_im,
                                                 self.bg_imgs.curr_stamp,
                                                 self.start_frame)
            
            #print 'time to perform background subtraction: '+str(time.time() - last_time)

            # process gui events
            if params.interactive:
                wx.Yield()
            if self.break_flag:
                break

            # find observations
            self.ellipses = ell.find_ellipses( self.dfore, self.isfore )

            #if params.DOBREAK:
            #    print 'Exiting at frame %d'%self.start_frame
            #    sys.exit(1)

            # process gui events
            if params.interactive:
                wx.Yield()
            if self.break_flag:
                break

            # match target identities to observations
            if len( self.ann_data ) > 1:
                flies = ell.find_flies( self.ann_data[-2],
                                        self.ann_data[-1],
                                        self.ellipses )
            elif len( self.ann_data ) == 1:
                flies = ell.find_flies( self.ann_data[0],
                                        self.ann_data[0],
                                        self.ellipses )
            else:
                flies = ell.TargetList()
                for i,obs in enumerate(self.ellipses):
                    if obs.isEmpty():
                        print 'empty observation'
                    else:
                        obs.identity = params.nids
                        flies.append(obs)
                        params.nids+=1

            # save to ann_data
            self.ann_data.append( flies )

            # fix any errors using hindsight
            self.hindsight.fixerrors()
            #print 'time to fix errors: '+str(time.time() - last_time)

            # write to file
            if (self.ann_data is not None) and \
                   (len(self.ann_data) > self.maxlookback+self.lastframewritten):
                self.lastframewritten = self.lastframewritten + 1
                self.ann_file.write_ellipses(self.ann_data[self.lastframewritten])
            #print 'In Track, after writing: self.start_frame = %d, params.start_frame = %d, len(ann_data) = %d'%(self.start_frame,params.start_frame,len(self.ann_data))

            # draw?
            if params.request_refresh or (params.do_refresh and ((self.start_frame % params.framesbetweenrefresh) == 0)):
                if params.interactive:
                    if self.start_frame:
                        self.ShowCurrentFrame()
                else:
                    sys.stderr.write("Frame %d / %d\n"%(self.start_frame,self.n_frames))
                params.request_refresh = False

            # process gui events
            if params.interactive:
                wx.Yield()
            if self.break_flag:
                break


        self.Finish()

    def Finish(self):

        # write the rest of the frames to file
        for ff in range(self.lastframewritten+1,len(self.ann_data)):
            self.ann_file.write_ellipses(self.ann_data[ff])
        self.lastframewritten = len(self.ann_data)-1
        # close the file
        self.ann_file.file.close()

    # enddef: Track()

    def StopThreads( self ):

        #wx.Yield()

        # stop algorithm
        self.break_flag = True

    def InitTrackingState(self):

	# initialize tracking state
        self.tracking = True
        self.ann_data = []
        params.nids = 0
        params.start_frame = 0
        self.start_frame = 0

    def DoAllPreprocessing(self):

        # estimate the background
        if (not self.IsBGModel()) or params.batch_autodetect_bg_model:
            print "Estimating background model"
            if params.interactive:
                self.bg_imgs.est_bg(self.frame)
            else:
                self.bg_imgs.est_bg()
        else:
            print "Not estimating background model"

        # detect arena if it has not been set yet
        if params.do_set_circular_arena and params.batch_autodetect_arena:
            print "Auto-detecting circular arena"
            setarena.doall(self.bg_imgs.center)
        else:
            print "Not detecting arena"

        self.bg_imgs.UpdateIsArena()

        # estimate the shape
        if params.batch_autodetect_shape:
            print "Estimating shape model"
            if params.interactive:
                ell.est_shape(self.bg_imgs,self.frame)
            else:
                ell.est_shape(self.bg_imgs)
        else:
            print "Not estimating shape model"

    def DoAll(self):

        if not params.interactive:
            self.RestoreStdio()

        sys.stderr.write("Performing preprocessing...\n")
	self.DoAllPreprocessing()

        sys.stderr.write("Done preprocessing, beginning tracking...\n")

        # begin tracking
        if params.interactive:
            self.UpdateToolBar('started')
        sys.stderr.write("Tracking...\n")
        self.Track()

        sys.stderr.write("Choosing Orientations...\n")
        # choose orientations
        self.choose_orientations = chooseorientations.ChooseOrientations(self.frame,self.ann_data,interactive=False)
        self.choose_orientations.ChooseOrientations()

        # save to a .mat file
        (basename,ext) = os.path.splitext(self.filename)
        if self.matfilename is None:
            savename = basename + '.mat'
        else:
            savename = self.matfilename
        sys.stderr.write("Saving to mat file "+savename+"...\n")

        self.ann_file.WriteMAT( savename, self.ann_data )
        sys.stderr.write("Done\n")
