# for some reason, tracking_settings does not work with wx2.8, only 2.6
#WXVER = '2.6'
#import wxversion
#if wxversion.checkInstalled(WXVER):
#  wxversion.select(WXVER)
USEGL = False
import wx
import params
import ellipsesk as ell
from wx import xrc
if USEGL:
  import motmot.wxglvideo.simple_overlay as wxvideo
else:
  import motmot.wxvideo.wxvideo as wxvideo
import motmot.wxvalidatedtext.wxvalidatedtext as wxvt # part of Motmot
import numpy as num
#from scipy.misc.pilutil import imresize
import imagesk
import copy
from matchidentities import cvpred

import pkg_resources # part of setuptools
RSRC_FILE = pkg_resources.resource_filename( __name__, "tracking_settings.xrc" )

SHOW_UNFILTERED_OBSERVATIONS = 0
SHOW_FILTERED_OBSERVATIONS = 1
SHOW_SMALL_OBSERVATIONS = 2
SHOW_LARGE_OBSERVATIONS = 3
SHOW_DELETED_OBSERVATIONS = 4
SHOW_SPLIT_OBSERVATIONS = 5
SHOW_MERGED_OBSERVATIONS = 6
SHOW_LOWERED_OBSERVATIONS = 7
SHOW_MAXJUMP = 8
SHOW_MOTIONMODEL = 9

class StoredObservations:

    def __init__(self,obs,frame,issmall=None,islarge=None,didlowerthresh=None,
                 didmerge=None,diddelete=None,didsplit=None):

        self.obs = obs
        self.frame = frame
        self.params = params.params.copy()
        self.issmall = issmall
        self.islarge = islarge
        self.didlowerthresh = didlowerthresh
        self.didmerge = didmerge
        self.diddelete = diddelete
        self.didsplit = didsplit

    def issame(self,frame):
        return ((self.params == params.params) and (self.frame == frame))

class TrackingSettings:

    def __init__(self,parent,bg_imgs,currframe):

        self.parent = parent
        # background model
        self.bg_imgs = bg_imgs
        self.show_frame = currframe
        self.bg_img_frame = -1
        # zoom mode
        self.zoomin = False
        self.zoomout = False
        self.info = False
        self.zoomaxes = [0,params.params.movie.get_width()-1,0,params.params.movie.get_height()-1]
        self.zoomfactor = 1

        self.shape_uptodate = False
        self.automatic_minshape = params.ShapeParams()
        self.automatic_maxshape = params.ShapeParams()
        self.automatic_meanshape = params.ShapeParams()

        rsrc = xrc.XmlResource( RSRC_FILE )
        self.frame = rsrc.LoadFrame(parent,"trackingframe")

        self.InitControlHandles()
        self.InitializeValues()
        self.BindCallbacks()
        self.OnResize()
        self.ShowImage()

    def InitControlHandles(self):

        #self.min_ntargets_input = self.control('min_ntarets')
        #self.max_ntargets_input = self.control('max_ntargets')
        self.automatic_shape_input = self.control('automatic_bounds')
        self.automatic_panel = self.control('automatic_panel')
        self.nstd_shape_input = self.control('nstd_shape')
        self.nframes_shape_input = self.control('nframes_shape')
        self.compute_shape_input = self.control('compute_shape')
        self.automatic_shape_text = self.control('automatic_shape_text')
        self.manual_shape_input = self.control('manual_bounds')
        self.manual_panel = self.control('manual_panel')
        self.min_area_input = self.control('min_area')
        self.max_area_input = self.control('max_area')
        self.min_major_input = self.control('min_major')
        self.max_major_input = self.control('max_major')
        self.min_minor_input = self.control('min_minor')
        self.max_minor_input = self.control('max_minor')
        self.min_ecc_input = self.control('min_ecc')
        self.max_ecc_input = self.control('max_ecc')
        self.angle_weight_input = self.control('angle_weight')
        self.max_jump_input = self.control('max_jump')
        self.center_dampen_input = self.control('center_dampen')
        self.angle_dampen_input = self.control('angle_dampen')
        self.max_area_delete_input = self.control('max_area_delete')
        self.min_area_ignore_input = self.control('min_area_ignore')
        self.max_penalty_merge_input = self.control('max_penalty_merge')
        self.lower_thresh_input = self.control('lower_thresh')
        self.done_input = self.control('done')
        self.img_panel = self.control('img_panel')
        self.frame_scrollbar = self.control('frame_scrollbar')
        self.frame_number_text = self.control('frame_number')
        self.show_img_input = self.control('show_img')
        #self.toolbar = self.frame.GetToolBar()
        self.toolbar = xrc.XRCCTRL(self.frame,'toolbar')
        self.zoomin_id = xrc.XRCID('zoomin')
        self.zoomout_id = xrc.XRCID('zoomout')
        self.info_id = xrc.XRCID('moreinfo')
        self.toolbar.AddSeparator()
        self.info_text = wx.TextCtrl(self.toolbar, -1, 'Observation Info', size=(300,20),style=wx.TE_READONLY|wx.TE_CENTRE)
        self.info_text.SetValue('Observation Info')
        #self.info_text.SetEditable(False)
        self.toolbar.AddControl(self.info_text)
        #self.zoomin_id = 1
        #self.zoomout_id = 2
        #self.info_id = 3
        self.hindsight_panel = self.control('hindsight_panel')
        self.splitdetection_panel = self.control('splitdetection_panel')
        self.mergeddetection_panel = self.control('mergeddetection_panel')
        self.spuriousdetection_panel = self.control('spuriousdetection_panel')
        self.lostdetection_panel = self.control('lostdetection_panel')
        self.do_fix_split_input = self.control('do_fix_split')
        self.do_fix_merged_input = self.control('do_fix_merged')
        self.do_fix_spurious_input = self.control('do_fix_spurious')
        self.do_fix_lost_input = self.control('do_fix_lost')
        self.splitdetection_length_input = self.control('splitdetection_length')
        self.splitdetection_cost_input = self.control('splitdetection_distance')
        self.mergeddetection_length_input = self.control('mergeddetection_length')
        self.mergeddetection_distance_input = self.control('mergeddetection_distance')
        self.spuriousdetection_length_input = self.control('spuriousdetection_length')
        self.lostdetection_length_input = self.control('lostdetection_length')

        box = wx.BoxSizer( wx.VERTICAL )
        self.img_panel.SetSizer( box )
        self.img_wind = wxvideo.DynamicImageCanvas( self.img_panel, -1 )
        self.img_wind.set_resize(True)
        box.Add( self.img_wind, 1, wx.EXPAND )
        self.img_panel.SetAutoLayout( True )
        self.img_panel.Layout()

    def InitializeValues(self):
        #self.min_ntargets_input.SetValue(str(params.params.min_ntargets))
        #self.max_ntargets_input.SetValue(str(params.params.max_ntargets))
        self.automatic_shape_input.SetValue(True)
        self.manual_shape_input.SetValue(False)
        self.nstd_shape_input.SetValue(str(params.params.n_std_thresh))
        self.nframes_shape_input.SetValue(str(params.params.n_frames_size))
        self.min_area_input.SetValue(str(params.params.minshape.area))
        self.min_major_input.SetValue(str(params.params.minshape.major))
        self.min_minor_input.SetValue(str(params.params.minshape.minor))
        self.min_ecc_input.SetValue(str(params.params.minshape.ecc))
        self.max_area_input.SetValue(str(params.params.maxshape.area))
        self.max_major_input.SetValue(str(params.params.maxshape.major))
        self.max_minor_input.SetValue(str(params.params.maxshape.minor))
        self.max_ecc_input.SetValue(str(params.params.maxshape.ecc))
        self.angle_dampen_input.SetValue(str(params.params.ang_dist_wt))
        self.max_area_delete_input.SetValue(str(params.params.maxareadelete))
        self.min_area_ignore_input.SetValue(str(params.params.minareaignore))
        self.max_penalty_merge_input.SetValue(str(params.params.maxpenaltymerge))
        self.max_jump_input.SetValue(str(params.params.max_jump))
        self.angle_weight_input.SetValue(str(params.params.ang_dist_wt))
        self.center_dampen_input.SetValue(str(params.params.dampen))
        self.angle_dampen_input.SetValue(str(params.params.angle_dampen))
        self.lower_thresh_input.SetValue(str(params.params.minbackthresh))
        self.frame_scrollbar.SetThumbPosition(self.show_frame)
        self.frame_scrollbar.SetScrollbar(self.show_frame,0,params.params.n_frames-1,30)
        self.img_chosen = SHOW_UNFILTERED_OBSERVATIONS
        self.show_img_input.SetSelection(self.img_chosen)
        self.toolbar.SetToggle(self.zoomin,self.zoomin_id)
        self.toolbar.SetToggle(self.zoomout,self.zoomout_id)
        self.toolbar.SetToggle(self.info,self.info_id)

        self.do_fix_split_input.SetValue(params.params.do_fix_split)
        self.do_fix_merged_input.SetValue(params.params.do_fix_merged)
        self.do_fix_spurious_input.SetValue(params.params.do_fix_spurious)
        self.do_fix_lost_input.SetValue(params.params.do_fix_lost)
        self.splitdetection_panel.Enable(params.params.do_fix_split)
        self.mergeddetection_panel.Enable(params.params.do_fix_merged)
        self.spuriousdetection_panel.Enable(params.params.do_fix_spurious)
        self.lostdetection_panel.Enable(params.params.do_fix_lost)
        
        self.splitdetection_length_input.SetValue(str(params.params.splitdetection_length))
        self.splitdetection_cost_input.SetValue('%.2f'%params.params.splitdetection_cost)
        self.mergeddetection_length_input.SetValue(str(params.params.mergeddetection_length))
        self.mergeddetection_distance_input.SetValue('%.2f'%params.params.mergeddetection_distance)
        self.spuriousdetection_length_input.SetValue(str(params.params.spuriousdetection_length))
        self.lostdetection_length_input.SetValue(str(params.params.lostdetection_length))

        
    def BindCallbacks(self):

        # number of targets
        #self.bindctrl(('min_ntargets','max_ntargets'),('int','int'),self.SetNTargets)

        # shape
        
        # automatic computation of bounds on shape
        self.frame.Bind(wx.EVT_RADIOBUTTON,self.OnAutomatic,self.automatic_shape_input)
        self.frame.Bind(wx.EVT_RADIOBUTTON,self.OnAutomatic,self.manual_shape_input)

        # automatic shape
        self.bindctrl(('nstd_shape','nframes_shape'),('float','float'),self.SetAutomatic)

        # compute shape now
        self.frame.Bind(wx.EVT_BUTTON,self.ComputeShapeNow,self.compute_shape_input)

        # manual shape
        self.bindctrl(('min_area','max_area','min_major','max_major',
                       'min_minor','max_minor','min_ecc','max_ecc'),
                      ('float','float','float','float',
                       'float','float','float','float'),
                      self.SetManual)

        # motion

        self.bindctrl(('angle_weight','max_jump','center_dampen','angle_dampen'),
                      ('float','float','float','float','float'),
                      self.SetMotion)

        # observation

        self.bindctrl(('max_area_delete','min_area_ignore','max_penalty_merge','lower_thresh'),
                      ('float','float','float','float'),self.SetObservation)


        # hindsight
        self.frame.Bind(wx.EVT_CHECKBOX,self.OnSplit,self.do_fix_split_input)
        self.frame.Bind(wx.EVT_CHECKBOX,self.OnMerged,self.do_fix_merged_input)
        self.frame.Bind(wx.EVT_CHECKBOX,self.OnSpurious,self.do_fix_spurious_input)
        self.frame.Bind(wx.EVT_CHECKBOX,self.OnLost,self.do_fix_lost_input)
        self.bindctrl(('splitdetection_length','splitdetection_distance'),
                      ('int','float'),self.SetSplit)
        self.bindctrl(('mergeddetection_length','mergeddetection_distance'),
                      ('int','float'),self.SetMerged)
        self.bindctrl(('spuriousdetection_length',),('int',),self.SetSpurious)
        self.bindctrl(('lostdetection_length',),('int',),self.SetLost)
        
        # frame scrollbar
        self.frame.Bind(wx.EVT_SCROLL,self.FrameScrollbarMoved,self.frame_scrollbar)

        # img choice 
	self.frame.Bind(wx.EVT_CHOICE,self.ImageChosen,self.show_img_input)

        # resize
        self.frame.Bind( wx.EVT_SIZE, self.OnResize )
        #self.frame.Bind( wx.EVT_MOVE, self.OnResizeBG )

        # toolbar
        self.frame.Bind(wx.EVT_TOOL, self.ZoominToggle, id=self.zoomin_id)
        self.frame.Bind(wx.EVT_TOOL, self.ZoomoutToggle, id=self.zoomout_id)
        self.frame.Bind(wx.EVT_TOOL, self.InfoToggle, id=self.info_id)

        # mouse click
        # get a pointer to the "trackset" child
        self.img_size = [params.params.movie.get_height(),params.params.movie.get_width()]
        img = num.zeros((self.img_size[0],self.img_size[1]),dtype=num.uint8)
        self.img_wind.update_image_and_drawings("trackset",
                                                img,
                                                format="MONO8")
        self.img_wind_child = self.img_wind.get_child_canvas("trackset")
        self.img_wind_child.Bind(wx.EVT_LEFT_DOWN,self.MouseClick)
        self.img_wind_child.Bind(wx.EVT_LEFT_DCLICK,self.MouseDoubleClick)

    def control(self,ctrlname):
        return xrc.XRCCTRL(self.frame,ctrlname)

    def bindctrl(self,ctrlnames,type,validatef):
        for i,v in enumerate(ctrlnames):
            if type[i] == 'int':
                wxvt.setup_validated_integer_callback(self.control(v),
                                                      xrc.XRCID(v),
                                                      validatef,
                                                      pending_color=params.params.wxvt_bg)
            else:
                wxvt.setup_validated_float_callback(self.control(v),
                                                    xrc.XRCID(v),
                                                    validatef,
                                                    pending_color=params.params.wxvt_bg)
    #def SetNTargets(self,evt):
    #    
    #    n1 = int(self.min_ntargets_input.GetValue())
    #    n2 = int(self.max_ntargets_input.GetValue())
    #
    #    if n1 <= n2:
    #        params.params.min_ntargets = n1
    #        params.params.max_ntargets = n2
    #    else:
    #        self.min_ntargets_input.SetValue(str(params.params.min_ntargets))
    #        self.max_ntargets_input.SetValue(str(params.params.max_ntargets))
    #
    #    self.ShowImage()

    def OnAutomatic(self,evt):

        isautomatic = self.automatic_shape_input.GetValue()
        if isautomatic == self.manual_shape_input.GetValue():
            print 'error: isautomatic == ismanual'
        if isautomatic:
            self.automatic_panel.Enable(True)
            self.manual_panel.Enable(False)
            if not self.shape_uptodate:
                self.automatic_shape_text.SetLabel('Bounds on shape not up to date.')
            else:
                self.automatic_shape_text.SetLabel('')
                params.params.minshape = self.automatic_minshape.copy()
                params.params.maxshape = self.automatic_maxshape.copy()
                params.params.meanshape = self.automatic_meanshape.copy()
                self.PrintShape()
                self.ShowImage()
        else:
            self.automatic_panel.Enable(False)
            self.manual_panel.Enable(True)
            self.automatic_shape_text.SetLabel('')

    def OnSplit(self,evt):

        params.params.do_fix_split = self.do_fix_split_input.GetValue()
        self.splitdetection_panel.Enable(params.params.do_fix_split)

    def OnMerged(self,evt):

        params.params.do_fix_merged = self.do_fix_merged_input.GetValue()
        self.mergeddetection_panel.Enable(params.params.do_fix_merged)

    def OnSpurious(self,evt):

        params.params.do_fix_spurious = self.do_fix_spurious_input.GetValue()
        self.spuriousdetection_panel.Enable(params.params.do_fix_spurious)

    def OnLost(self,evt):

        params.params.do_fix_lost = self.do_fix_lost_input.GetValue()
        self.lostdetection_panel.Enable(params.params.do_fix_lost)

    def SetSplit(self,evt):

        params.params.splitdetection_length = int(self.splitdetection_length_input.GetValue())
        params.params.splitdetection_cost = float(self.splitdetection_cost_input.GetValue())

    def SetMerged(self,evt):

        params.params.mergeddetection_length = int(self.mergeddetection_length_input.GetValue())
        params.params.mergeddetection_distance = float(self.mergeddetection_distance_input.GetValue())

    def SetSpurious(self,evt):

        params.params.spuriousdetection_length = int(self.spuriousdetection_length_input.GetValue())

    def SetLost(self,evt):

        params.params.lostdetection_length = int(self.lostdetection_length_input.GetValue())

    def SetAutomatic(self,evt):
        nstd = float(self.nstd_shape_input.GetValue())
        nframes = int(self.nframes_shape_input.GetValue())

        # check to make sure valid
        if nstd <= 0:
            self.nstd_shape_input.SetValue(str(params.params.n_std_thresh))
            nstd = params.params.n_std_thresh

        if nframes <= 0:
            self.nframes_shape_input.SetValue(str(params.params.n_frames_size))
            nframes = params.params.n_frames_size
        
        if (nstd != params.params.n_std_thresh) or (nframes != params.params.n_frames_size):
            params.params.n_std_thresh = nstd
            params.params.n_frames_size = nframes
            self.shape_uptodate = False
            self.automatic_shape_text.SetLabel('Bounds on shape not up to date.')

    def ComputeShapeNow(self,evt):

        # estimate shape now
        wx.BeginBusyCursor()
        wx.Yield()
        ell.est_shape( self.bg_imgs )
        wx.EndBusyCursor()
        
        # copy to temporary variable
        self.automatic_minshape = params.params.minshape.copy()
        self.automatic_maxshape = params.params.maxshape.copy()
        self.automatic_meanshape = params.params.meanshape.copy()

        # show shape
        self.PrintShape()

        # set up to date
        self.shape_uptodate = True
        self.automatic_shape_text.SetLabel('')

        self.ShowImage()

    def PrintShape(self):

        self.min_area_input.SetValue(str(params.params.minshape.area))
        self.min_major_input.SetValue(str(params.params.minshape.major))
        self.min_minor_input.SetValue(str(params.params.minshape.minor))
        self.min_ecc_input.SetValue(str(params.params.minshape.ecc))
        self.max_area_input.SetValue(str(params.params.maxshape.area))
        self.max_major_input.SetValue(str(params.params.maxshape.major))
        self.max_minor_input.SetValue(str(params.params.maxshape.minor))
        self.max_ecc_input.SetValue(str(params.params.maxshape.ecc))

    def SetManual(self,evt):

        minarea = float(self.min_area_input.GetValue())
        maxarea = float(self.max_area_input.GetValue())
        if (minarea < maxarea) and (maxarea > 0):
            params.params.minshape.area = minarea
            params.params.maxshape.area = maxarea
        else:
            self.min_area_input.SetValue(str(params.params.minshape.area))
            self.max_area_input.SetValue(str(params.params.maxshape.area))
        minmajor = float(self.min_major_input.GetValue())
        maxmajor = float(self.max_major_input.GetValue())
        if (minmajor < maxmajor) and (maxmajor > 0):
            params.params.minshape.major = minmajor
            params.params.maxshape.major = maxmajor
        else:
            self.min_major_input.SetValue(str(params.params.minshape.major))
            self.max_major_input.SetValue(str(params.params.maxshape.major))
        minminor = float(self.min_minor_input.GetValue())
        maxminor = float(self.max_minor_input.GetValue())
        if (minminor < maxminor) and (maxminor > 0):
            params.params.minshape.minor = minminor
            params.params.maxshape.minor = maxminor
        else:
            self.min_minor_input.SetValue(str(params.params.minshape.minor))
            self.max_minor_input.SetValue(str(params.params.maxshape.minor))
        minecc = float(self.min_ecc_input.GetValue())
        maxecc = float(self.max_ecc_input.GetValue())
        if (minecc < maxecc) and (maxecc > 0) and (minecc < 1):
            params.params.minshape.ecc = minecc
            params.params.maxshape.ecc = maxecc
        else:
            self.min_ecc_input.SetValue(str(params.params.minshape.ecc))
            self.max_ecc_input.SetValue(str(params.params.maxshape.ecc))

        params.params.meanshape = params.averageshape(params.params.minshape,
                                                      params.params.maxshape)

        self.ShowImage()

    def SetMotion(self,evt):
        
        angle_weight = float(self.angle_weight_input.GetValue())
        if angle_weight >= 0:
            params.params.ang_dist_wt = angle_weight
        else:
            self.angle_weight_input.SetValue(str(params.params.ang_dist_wt))

        max_jump = float(self.max_jump_input.GetValue())
        if max_jump > 0:
            params.params.max_jump = max_jump
        else:
            self.max_jump_input.SetValue(str(params.params.max_jump))

        center_dampen = float(self.center_dampen_input.GetValue())
        if center_dampen >= 0 and center_dampen <= 1:
            params.params.dampen = center_dampen
        else:
            self.center_dampen_input.SetValue(str(params.params.dampen))

        angle_dampen = float(self.angle_dampen_input.GetValue())
        if angle_dampen >= 0 and angle_dampen <= 1:
            params.params.angle_dampen = angle_dampen
        else:
            self.angle_dampen_input.SetValue(str(params.params.angle_dampen))

        self.ShowImage()

    def SetObservation(self,evt):

        params.params.maxareadelete = float(self.max_area_delete_input.GetValue())
        params.params.minareaignore = float(self.min_area_ignore_input.GetValue())
        params.params.maxpenaltymerge = float(self.max_penalty_merge_input.GetValue())
        minbackthresh = float(self.lower_thresh_input.GetValue())
        if minbackthresh > 0 and minbackthresh <= 1:
            params.params.minbackthresh = minbackthresh
        else:
            self.lower_thresh_input.GetValue(params.params.minbackthresh)

        self.ShowImage()

    def FrameScrollbarMoved(self,evt):
	
	# get the value on the scrollbar now
	self.show_frame = self.frame_scrollbar.GetThumbPosition()

        if hasattr(self,'obs_filtered'):
            self.obs_prev = self.obs_filtered

	# show the image and update text
	self.ShowImage()

    def ShowImage(self):

        im, stamp = params.params.movie.get_frame( int(self.show_frame) )
        windowsize = [self.img_panel.GetRect().GetHeight(),self.img_panel.GetRect().GetWidth()]

        self.GetBgImage()
            
        if self.img_chosen == SHOW_UNFILTERED_OBSERVATIONS:
            obs_unfiltered = self.GetObsUnfiltered()
            plot_linesegs = ell.draw_ellipses(obs_unfiltered)
            
        elif self.img_chosen == SHOW_FILTERED_OBSERVATIONS:
            obs_filtered = self.GetObsFiltered()
            plot_linesegs = ell.draw_ellipses(obs_filtered)
            
        elif self.img_chosen == SHOW_SMALL_OBSERVATIONS:
            obs_unfiltered = self.GetObsUnfiltered()
            obs_small = []
            for obs in obs_unfiltered:
                if obs.area < params.params.minshape.area:
                    obs_small.append(obs)
            plot_linesegs = ell.draw_ellipses(obs_small)
            
        elif self.img_chosen == SHOW_LARGE_OBSERVATIONS:
            obs_unfiltered = self.GetObsUnfiltered()
            obs_large = []
            for obs in obs_unfiltered:
                if obs.area > params.params.maxshape.area:
                    obs_large.append(obs)
            plot_linesegs = ell.draw_ellipses(obs_large)

        elif self.img_chosen == SHOW_DELETED_OBSERVATIONS:
            (obs_unfiltered,diddelete) = self.GetObsUnfiltered('diddelete')
            obs_deleted = []
            for (i,v) in enumerate(diddelete):
                if v: obs_deleted.append(obs_unfiltered[i])
            plot_linesegs = ell.draw_ellipses(obs_deleted)

        elif self.img_chosen == SHOW_SPLIT_OBSERVATIONS:
            (obs_unfiltered,didsplit) = self.GetObsUnfiltered('didsplit')
            colors = []
            obs_split = []
            for (i,v) in enumerate(didsplit):
                if len(v)>1:
                    for j in v:
                        colors.append(params.params.colors[i % len(params.params.colors)])
                        obs_split.append(obs_unfiltered[j])
            plot_linesegs = ell.draw_ellipses(obs_split)

        elif self.img_chosen == SHOW_MERGED_OBSERVATIONS:
            (obs_unfiltered,didmerge) = self.GetObsUnfiltered('didmerge')
            colors = []
            obs_merge = []
            for (i,v) in enumerate(didmerge):
                if len(v)>1:
                  obs_merge.append(obs_unfiltered[i])
            plot_linesegs = ell.draw_ellipses(obs_merge)

        elif self.img_chosen == SHOW_LOWERED_OBSERVATIONS:
            (obs_unfiltered,didlowerthresh) = self.GetObsUnfiltered('didlowerthresh')
            obs_lowered = []
            for (i,v) in enumerate(didlowerthresh):
                if v: obs_lowered.append(obs_unfiltered[i])
            plot_linesegs = ell.draw_ellipses(obs_lowered)

        # MAXJUMP
        elif self.img_chosen == SHOW_MAXJUMP:

            # either grab or compute observations
            obs_filtered = self.GetObsFiltered()
            plot_linesegs = ell.draw_ellipses(obs_filtered)

            # draw circles
            for i,obs in enumerate(obs_filtered):
                plot_new_stuff = imagesk.draw_circle(obs.center.x,
                                                     obs.center.y,params.params.max_jump,
                                                     params.params.colors[i%len(params.params.colors)])
                plot_linesegs.extend(plot_new_stuff)
                

        elif self.img_chosen == SHOW_MOTIONMODEL:

            (target_prev,target_curr,target_pred) = self.GetTargetMotion()

            # show the next frame, if there is one
            nextframe = num.minimum(int(self.show_frame+1),params.params.n_frames-1)
            im, stamp = params.params.movie.get_frame(nextframe)

            # draw previous positions
            plot_linesegs = ell.draw_ellipses(target_prev)

            # draw current positions
            plot_new_stuff = ell.draw_ellipses(target_curr,colors=[[255,255,0]])
            plot_linesegs.extend(plot_new_stuff)

            # draw predicted positions
            plot_new_stuff = ell.draw_ellipses(target_pred)
            plot_linesegs.extend(plot_new_stuff)

            scaleunit = params.params.max_jump / params.params.DRAW_MOTION_SCALE
            parcolor = [0,255,0]
            perpcolor = [0,255,0]
            anglecolor = [0,0,255]

            for i in target_pred.iterkeys():
                # compute direction of motion
                vx = target_pred[i].center.x - target_curr[i].center.x
                vy = target_pred[i].center.y - target_curr[i].center.y
                thetamotion = num.arctan2(vy,vx)

                # angle
                theta = target_pred[i].angle
                
                # we want to choose to add pi if that makes difference from motion direction smaller
                dtheta = abs( ((theta - thetamotion + num.pi) % (2.*num.pi)) - num.pi )
                if dtheta > (num.pi/2.):
                    theta += num.pi

                # compute end points of parallel motion
                x0 = target_pred[i].center.x + scaleunit*num.cos(theta)
                x1 = target_pred[i].center.x - scaleunit*num.cos(theta)
                y0 = target_pred[i].center.y + scaleunit*num.sin(theta)
                y1 = target_pred[i].center.y - scaleunit*num.sin(theta)

                # add parallel motion annotation line
                plot_new_stuff = imagesk.draw_line(x0+1,y0+1,x1+1,y1+1,parcolor)
                plot_linesegs.extend(plot_new_stuff)

                # compute end points of perpendicular motion
                x0 = target_pred[i].center.x + scaleunit*num.cos(num.pi/2.+theta)
                x1 = target_pred[i].center.x - scaleunit*num.cos(num.pi/2.+theta)
                y0 = target_pred[i].center.y + scaleunit*num.sin(num.pi/2.+theta)
                y1 = target_pred[i].center.y - scaleunit*num.sin(num.pi/2.+theta)
                
                # add perpendicular motion annotation line
                plot_new_stuff = imagesk.draw_line(x0+1,y0+1,x1+1,y1+1,perpcolor)
                plot_linesegs.extend(plot_new_stuff)
                

                # compute end points of angular motion
                if params.params.ang_dist_wt > 0:
                    
                    dtheta = scaleunit*num.sqrt(1./params.params.ang_dist_wt)
                    if dtheta >= (num.pi/2.):
                        print 'dtheta is more than pi/2'
                    else:
                        theta0 = theta - dtheta
                        theta1 = theta + dtheta
                        
                        # draw arc
                        plot_new_stuff = imagesk.draw_arc(target_pred[i].center.x,
                                                          target_pred[i].center.y,
                                                          scaleunit/2,
                                                          theta0,theta1,
                                                          anglecolor)
                        plot_linesegs.extend(plot_new_stuff)
                   
        im = imagesk.double2mono8(im,donormalize=False)
        linesegs,im = imagesk.zoom_linesegs_and_image(plot_linesegs,im,self.zoomaxes)
        (linesegs,linecolors) = imagesk.separate_linesegs_colors(linesegs)
        
        self.img_wind.update_image_and_drawings('trackset',
                                                im,
                                                format='MONO8',
                                                linesegs=linesegs,
                                                lineseg_colors=linecolors
                                                )

        self.img_wind.Refresh(eraseBackground=False)

        self.frame_number_text.SetLabel('Frame %d'%self.show_frame)

    def ImageChosen(self,evt):

        self.img_chosen = self.show_img_input.GetSelection()
        self.ShowImage()

    def OnResize(self,evt=None):

        if evt is not None: evt.Skip()
        self.frame.Layout()
        try:
            #self.ShowImage()
            self.frame_scrollbar.SetMinSize( wx.Size(
                self.img_panel.GetRect().GetWidth(),
                self.frame_scrollbar.GetRect().GetHeight() ) )
        except AttributeError: pass # during initialization

    def ZoominToggle(self,evt):

        self.zoomin = self.zoomin == False
        if self.zoomin == True and self.zoomout == True:
            self.toolbar.ToggleTool(self.zoomout_id,False)
            self.zoomout = False
        if self.zoomin == True and self.info == True:
            self.toolbar.ToggleTool(self.info_id,False)
            self.info = False

    def ZoomoutToggle(self,evt):

        self.zoomout = self.zoomout == False
        if self.zoomin == True and self.zoomout == True:
            self.toolbar.ToggleTool(self.zoomin_id,False)
            self.zoomin = False
        if self.zoomout == True and self.info == True:
            self.toolbar.ToggleTool(self.info_id,False)
            self.info = False

    def InfoToggle(self,evt):

        self.info = self.info == False
        if self.info == True and self.zoomin == True:
            self.toolbar.ToggleTool(self.zoomin_id,False)
            self.zoomin = False
        if self.zoomout == True and self.info == True:
            self.toolbar.ToggleTool(self.zoomout_id,False)
            self.zoomout = False

    def MouseDoubleClick(self,evt):

        if self.zoomout == True:
            self.zoomfactor = 1
            self.SetZoomAxes()

    def MouseClick(self,evt):

        # get the clicked position
        if USEGL:
          windowheight = self.img_wind_child.GetRect().GetHeight()
          windowwidth = self.img_wind_child.GetRect().GetWidth()
          x = evt.GetX() * self.img_size[1] / windowwidth
          y = self.img_size[0] - evt.GetY() * self.img_size[0] / windowheight
        else:
          resize = self.img_wind.get_resize()
          x = evt.GetX() / resize
          h = self.zoomaxes[3] - self.zoomaxes[2] + 1
          y = h - evt.GetY() / resize
          
        x += self.zoomaxes[0]
        y += self.zoomaxes[2]

        if (x > self.img_size[1]) or (y > self.img_size[0]):
            return
        if self.zoomin:
            self.ZoomIn(x,y)
        elif self.zoomout:
            self.ZoomOut()
        elif self.info:
            self.ShowInfo(x,y)
        
    def ZoomIn(self,x,y):

        self.zoomfactor *= 1.5
        self.zoompoint = [x,y]
        self.SetZoomAxes()

    def SetZoomAxes(self):

        x = self.zoompoint[0]
        y = self.zoompoint[1]
        
        W = params.params.movie.get_width()
        H = params.params.movie.get_height()
        h = H/self.zoomfactor
        w = W/self.zoomfactor
        x1 = x-w/2
        x2 = x+w/2
        y1 = y-h/2
        y2 = y+h/2

        if x1 < 0:
            x2 -= x1
            x1 = 0
        elif x2 > W-1:
            x1 -= (x2 - W + 1)
            x2 = W-1
        if y1 < 0:
            y2 -= y1
            y1 = 0
        elif y2 > H-1:
            y1 -= (y2 - H + 1)
            y2 = H-1
        x1 = num.maximum(int(x1),0)
        x2 = num.minimum(int(x2),W-1)
        y1 = num.maximum(int(y1),0)
        y2 = num.minimum(int(y2),H-1)

        self.zoomaxes = [x1,x2,y1,y2]
        self.ShowImage()

    def ZoomOut(self):

        if self.zoomfactor <= 1:
            return

        self.zoomfactor /= 1.5
        self.SetZoomAxes()

    def ShowInfo(self,x,y):

        # grab targets
        if (self.img_chosen == SHOW_MAXJUMP) or (self.img_chosen == SHOW_FILTERED_OBSERVATIONS):
            obs = self.obs_filtered.obs
        else:
            obs = self.obs_unfiltered.obs

        # determine closest target
        mind = num.inf
        for i,v in enumerate(obs):
            d = (v.center.x-x)**2 + (v.center.y-y)**2
            if d <= mind:
                mini = i
                mind = d

        maxdshowinfo = (num.maximum(self.zoomaxes[1]-self.zoomaxes[0],
                                   self.zoomaxes[2]-self.zoomaxes[1])/params.params.MAXDSHOWINFO)**2

        if mind < maxdshowinfo:
            self.ShowObsInfo(obs[mini])
        
    def ShowObsInfo(self,ellipse):
        self.info_text.SetValue('area=%.2f, maj=%.2f, min=%.2f, ecc=%.2f'%(ellipse.area,ellipse.major,ellipse.minor,ellipse.minor/ellipse.major))

    def GetObsFiltered(self):
        if not(hasattr(self,'obs_filtered') and self.obs_filtered.issame(self.show_frame)):
            wx.BeginBusyCursor()
            wx.Yield()
            obs_filtered = ell.find_ellipses(self.bg_imgs.dfore,self.bg_imgs.bw,True)
            wx.EndBusyCursor()
            self.obs_filtered = StoredObservations(obs_filtered,self.show_frame)

        return self.obs_filtered.obs
    def GetObsUnfiltered(self,*args):

        # if we are only interested in the unfiltered observation
        if len(args) == 0:
            # if it has not yet been computed for this frame, compute
            if not(hasattr(self,'obs_unfiltered') and self.obs_unfiltered.issame(self.show_frame)):
                wx.BeginBusyCursor()
                wx.YieldIfNeeded()
                obs_unfiltered = ell.find_ellipses(self.bg_imgs.dfore,self.bg_imgs.bw,False)
                wx.EndBusyCursor()
                self.obs_unfiltered = StoredObservations(obs_unfiltered,self.show_frame)
            return self.obs_unfiltered.obs

        # otherwise, we need to check for attributes as well
        mustcompute = False
        if hasattr(self,'obs_unfiltered') and self.obs_unfiltered.issame(self.show_frame):
            for arg in args:
               if self.obs_unfiltered.__dict__[arg] is None:
                   mustcompute = True
                   break
        else:
            mustcompute = True
            
        # compute if necessary
        if mustcompute:
            wx.BeginBusyCursor()
            wx.YieldIfNeeded()
            (obs_unfiltered,issmall,islarge,didlowerthresh,didmerge,diddelete,didsplit) = \
                                         ell.find_ellipses_display(self.bg_imgs.dfore,self.bg_imgs.bw)
            wx.EndBusyCursor()
            self.obs_unfiltered = StoredObservations(obs_unfiltered,self.show_frame,
                                                     issmall,islarge,didlowerthresh,
                                                     didmerge,diddelete,didsplit)

        # create return list
        ret = (self.obs_unfiltered.obs,)
        for arg in args:
            ret += (self.obs_unfiltered.__dict__[arg],)

        return ret

    def GetBgImage(self):
        
        if not (self.bg_img_frame == self.show_frame):
            (self.bg_imgs.dfore,self.bg_imgs.bw) = self.bg_imgs.sub_bg(self.show_frame)
            self.bg_img_frame = self.show_frame
        
    def GetObsPrev(self):
        if not(hasattr(self,'obs_prev') and self.obs_filtered.issame(self.show_frame-1)):
            wx.BeginBusyCursor()
            wx.Yield()
            prevframe = num.maximum(0,self.show_frame-1)
            (dfore,bw) = self.bg_imgs.sub_bg(prevframe)
            obs_filtered = ell.find_ellipses(dfore,bw,True)
            wx.EndBusyCursor()
            self.obs_prev = StoredObservations(obs_filtered,self.show_frame)
        return self.obs_prev.obs

    def GetTargetMotion(self):

        # if it has already been computed, then just grab and return
        if hasattr(self,'targetmotion') and self.targetmotion[3] == self.show_frame \
           and self.targetmotion[4] == params.params:
            return self.targetmotion[0:3]
        
        # get current positions
        obs_curr = self.GetObsFiltered()
        # get previous positions
        obs_prev = self.GetObsPrev()
        # give identities to previous positions
        target_prev = ell.TargetList()
        for i,obs in enumerate(obs_prev):
            obs.identity = i
            target_prev.append(obs)
        # match previous and current targets, no velocity
        oldnids = params.params.nids
        target_curr = ell.find_flies(target_prev,target_prev,obs_curr)
        # don't actually assign new identities
        params.params.nids = oldnids
        # delete targets that aren't in both frames
        keyscurr = set(target_curr.keys())
        keysprev = set(target_prev.keys())
        keysremove = keyscurr - keysprev
        for i in keysremove:
            tmp = target_curr.pop(i)
        keysremove = keysprev - keyscurr
        for i in keysremove:
            tmp = target_prev.pop(i)
        
        # compute predicted positions
        target_pred = cvpred(target_prev,target_curr)
        # store
        self.targetmotion = (target_prev,target_curr,target_pred,self.show_frame,params.params.copy())
        # return
        return self.targetmotion[0:3]