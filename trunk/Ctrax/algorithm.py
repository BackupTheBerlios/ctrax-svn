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

import motmot.wxvalidatedtext as wxvt # part of Motmot

from version import DEBUG
import annfiles as annot
import ellipsesk as ell
import settings
import hindsight
import sys

from params import params

import pkg_resources # part of setuptools
SETTINGS_RSRC_FILE = pkg_resources.resource_filename( __name__, "track_settings.xrc" )

# CtraxApp.Track #################################################
class CtraxAlgorithm (settings.AppWithSettings):
    """Cannot be used alone -- this class only exists
    to keep algorithm code together in one file."""
    def Track( self ):
        """Run the m-tracker."""

        ## initialization ##

        if DEBUG: print "Tracking from frame %d..."%self.start_frame
        if DEBUG: print "Last frame tracked = %d"%self.ann_file.lastframetracked


        # maximum number of frames we will look back to fix errors
        self.maxlookback = max(params.lostdetection_length,
                               params.spuriousdetection_length,
                               params.mergeddetection_length,
                               params.splitdetection_length)

        if params.interactive:
            wx.Yield()

        # initialize hindsight data structures

        if DEBUG: print "Initializing hindsight structures"
        self.hindsight = hindsight.Hindsight(self.ann_file,self.bg_imgs)
        self.hindsight.initialize_milestones()
        if DEBUG: print "Done initializing hindsight structures"

        self.break_flag = False

        if DEBUG: print "Initializing buffer for tracking"
        self.ann_file.InitializeBufferForTracking(self.start_frame)

        for self.start_frame in range(self.start_frame,self.n_frames):

            if DEBUG: print "Tracking frame %d / %d"%(self.start_frame,self.n_frames-1)
        
            #if DEBUG:
            #    break

            if self.break_flag:
                break

            last_time = time.time()

            # perform background subtraction
            try:
                (self.dfore,self.isfore) = self.bg_imgs.sub_bg( self.start_frame )
            except:
                # catch all error types here, and just break out of loop
                break

            # write to sbfmf
            if self.dowritesbfmf:
                self.movie.writesbfmf_writeframe(self.isfore,
                                                 self.bg_imgs.curr_im,
                                                 self.bg_imgs.curr_stamp,
                                                 self.start_frame)
            
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
            if len( self.ann_file ) > 1:
                flies = ell.find_flies( self.ann_file[-2],
                                        self.ann_file[-1],
                                        self.ellipses,
                                        self.ann_file)
            elif len( self.ann_file ) == 1:
                flies = ell.find_flies( self.ann_file[0],
                                        self.ann_file[0],
                                        self.ellipses,
                                        self.ann_file)
            else:
                flies = ell.TargetList()
                for i,obs in enumerate(self.ellipses):
                    if obs.isEmpty():
                        if DEBUG: print 'empty observation'
                    else:
                        newid = self.ann_file.GetNewId()
                        obs.identity = newid
                        flies.append(obs)

            if DEBUG: print "Done with frame %d, appending to ann_file"%self.start_frame

            # save to ann_data
            self.ann_file.append( flies )

            if DEBUG: print "Added to ann_file, now running fixerrors"

            # fix any errors using hindsight
            self.hindsight.fixerrors()
            #print 'time to fix errors: '+str(time.time() - last_time)

            

            # draw?
            if params.request_refresh or (params.do_refresh and ((self.start_frame % params.framesbetweenrefresh) == 0)):
                if params.interactive:
                    if self.start_frame:
                        self.ShowCurrentFrame()
                else:
                    print "Frame %d / %d\n"%(self.start_frame,self.n_frames)
                params.request_refresh = False

            # process gui events
            if params.interactive:
                wx.Yield()
            if self.break_flag:
                break


        self.Finish()

    def Finish(self):

        # write the rest of the frames to file
        self.ann_file.finish_writing()

    # enddef: Track()

    def StopThreads( self ):

        #wx.Yield()

        # stop algorithm
        self.break_flag = True

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

        print "Performing preprocessing...\n"
	self.DoAllPreprocessing()

        print "Done preprocessing, beginning tracking...\n"

        # begin tracking
        if params.interactive:
            self.UpdateToolBar('started')

        # write sbfmf header
        if self.dowritesbfmf:
            # open an sbfmf file if necessary
            self.movie.writesbfmf_start(self.bg_imgs,
                                        self.writesbfmf_filename)

        print "Tracking..."
        try:
            self.Track()
        except:
            print "Error during Track"
            raise
        print "Done tracking"

        # write the sbfmf index and close the sbfmf file
        if self.dowritesbfmf and self.movie.writesbfmf_isopen():
            self.movie.writesbfmf_close(self.start_frame)
        
        print "Choosing Orientations..."
        # choose orientations
        self.choose_orientations = chooseorientations.ChooseOrientations(self.frame,interactive=False)
        self.ann_file = self.choose_orientations.ChooseOrientations(self.ann_file)

        # save to a .mat file
        (basename,ext) = os.path.splitext(self.filename)
        if not(hasattr(self,'matfilename')) or self.matfilename is None:
            savename = basename + '.mat'
        else:
            savename = self.matfilename
        print "Saving to mat file "+savename+"...\n"

        self.ann_file.WriteMAT( savename )
        print "Done\n"

