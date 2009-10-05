# Ctrax_settings.py
# KMB 09/12/07

import os
import sys
import wx
from wx import xrc
from params import params
import tracking_settings
import chooseorientations
import numpy as num

import motmot.wxglvideo.simple_overlay as wxvideo
import motmot.wxvalidatedtext.wxvalidatedtext as wxvt

from version import __version__, DEBUG

import pkg_resources # part of setuptools
RSRC_FILE = pkg_resources.resource_filename( __name__, "Ctrax.xrc" )

# these may import 'const'
import algorithm as alg
import bg
import draw
import ellipsesk as ell

ID_PRINTINTERVALS = 837487
ID_PRINTBUFFEREL = 837488
ID_PRINTFRAME = 837489

class AppWithSettings( wx.App ):
    """Cannot be used alone -- this class only exists
    to keep settings GUI code together in one file."""
    def InitGUI( self ):
        """Load XML resources, create handles, bind callbacks."""
        rsrc = xrc.XmlResource( RSRC_FILE )

        self.frame = rsrc.LoadFrame( None, "frame_Ctrax" )
        self.frame.SetIcon(wx.Icon(pkg_resources.resource_filename( __name__,'Ctraxicon.ico'), wx.BITMAP_TYPE_ICO))

        # make references to useful objects
        self.menu = self.frame.GetMenuBar()
        if DEBUG: 
            self.debugmenu = wx.Menu()
            self.debugmenu.Append(ID_PRINTINTERVALS,"Print Intervals",
                                  "Print frame intervals tracked, buffered, and written.")
            self.debugmenu.Append(ID_PRINTBUFFEREL,"Print Buffer...",
                                  "Print ellipses stored at specified element of buffer.")
            self.debugmenu.Append(ID_PRINTFRAME,"Print Frame...",
                                  "Pring ellipses for specified frame.")
            self.menu.Append( self.debugmenu, "DEBUG" )
            
        self.status = xrc.XRCCTRL( self.frame, "bar_status" )
        self.framenumber_text = xrc.XRCCTRL( self.frame, "text_framenumber" )
        self.num_flies_text = xrc.XRCCTRL( self.frame, "text_num_flies" )
        self.rate_text = xrc.XRCCTRL( self.frame, "text_refresh_rate")
        self.slider = xrc.XRCCTRL( self.frame, "slider_frame" )
        self.slider.Enable( False )

        # make image window
        self.img_panel = xrc.XRCCTRL( self.frame, "panel_img" )
        box = wx.BoxSizer( wx.VERTICAL )
        self.img_panel.SetSizer( box )
        self.img_wind = wxvideo.DynamicImageCanvas( self.img_panel, -1 )
        box.Add( self.img_wind, 1, wx.EXPAND )
        self.img_panel.SetAutoLayout( True )
        self.img_panel.Layout()

        self.zoommode = False

        # set up tools
        self.toolbar = xrc.XRCCTRL(self.frame,'toolbar')

        # other tools
        self.zoom_id = xrc.XRCID('zoom')
        self.play_id = xrc.XRCID('play')
        self.stop_id = xrc.XRCID('stop')
        self.speedup_id = xrc.XRCID('speed_up')
        self.slowdown_id = xrc.XRCID('slow_down')
        self.refresh_id = xrc.XRCID('refresh')

        # set up appearances for toolbar
        self.toolbar.SetToggle(self.zoommode,self.zoom_id)
        self.stop_tracking_tooltip = 'Stop Tracking'
        self.stop_playback_tooltip = 'Stop Video Playback'
        self.start_playback_tooltip = 'Start Video Playback'
        self.speedup_tracking_tooltip = 'Increase Refresh Rate'
        self.speedup_playback_tooltip = 'Increase Playback Speed'
        self.slowdown_tracking_tooltip = 'Decrease Refresh Rate'
        self.slowdown_playback_tooltip = 'Decrease Playback Speed'
        self.play_speed = 1.0

        self.UpdateToolBar('stopped')

        # bind callbacks
        # file menu
        self.frame.Bind( wx.EVT_MENU, self.OnOpen, id=xrc.XRCID("menu_file_open") )
        self.frame.Bind( wx.EVT_MENU, self.OnLoadSettings, id=xrc.XRCID("menu_load_settings") )
        self.frame.Bind( wx.EVT_MENU, self.OnBatch, id=xrc.XRCID("menu_file_batch") )
        self.frame.Bind( wx.EVT_MENU, self.OnSave, id=xrc.XRCID("menu_file_export") )
        self.frame.Bind( wx.EVT_MENU, self.OnSaveTimestamps, id=xrc.XRCID("menu_file_write_timestamps") )
        self.frame.Bind( wx.EVT_MENU, self.OnSaveAvi, id=xrc.XRCID("menu_file_save_avi") )
        self.frame.Bind( wx.EVT_MENU, self.OnQuit, id=xrc.XRCID("menu_file_quit") )
        # help menu
        self.frame.Bind( wx.EVT_MENU, self.OnHelp, id=xrc.XRCID("menu_help_help") )
        self.frame.Bind( wx.EVT_MENU, self.OnAbout, id=xrc.XRCID("menu_help_about") )
        # track menu
        self.frame.Bind( wx.EVT_MENU, self.OnStartTrackingMenu, id=xrc.XRCID("menu_track_start") )
        self.frame.Bind( wx.EVT_MENU, self.OnStartTrackingMenu, id=xrc.XRCID("menu_track_resume") )
        self.frame.Bind( wx.EVT_MENU, self.OnStartTrackingMenu, id=xrc.XRCID("menu_track_resume_here") )
        self.frame.Bind( wx.EVT_MENU, self.OnWriteSBFMF, id=xrc.XRCID("menu_track_writesbfmf") )
        self.frame.Bind( wx.EVT_MENU, self.OnComputeBg, id=xrc.XRCID("menu_compute_background") )
        self.frame.Bind( wx.EVT_MENU, self.OnComputeShape, id=xrc.XRCID("menu_compute_shape") )
        # settings menu
        self.frame.Bind( wx.EVT_MENU, self.OnSettingsBGModel, id=xrc.XRCID("menu_settings_bg_model") )
        self.frame.Bind( wx.EVT_MENU, self.OnSettingsBG, id=xrc.XRCID("menu_settings_bg") )
        self.frame.Bind( wx.EVT_MENU, self.OnSettingsTracking, id=xrc.XRCID("menu_settings_tracking") )
        self.frame.Bind( wx.EVT_MENU, self.OnChooseOrientations,id=xrc.XRCID("menu_choose_orientations"))
        #self.frame.Bind( wx.EVT_MENU, self.OnEllipseSize, id=xrc.XRCID("menu_settings_ellipses") )
        #self.frame.Bind( wx.EVT_MENU, self.OnCheckColorblind, id=xrc.XRCID("menu_settings_use_colorblind") )
        self.frame.Bind( wx.EVT_MENU, self.OnCheckShowAnn, id=xrc.XRCID("menu_playback_show_ann") )
        self.frame.Bind( wx.EVT_MENU, self.OnCheckRefresh, id=xrc.XRCID("menu_do_refresh") )
        self.frame.Bind( wx.EVT_MENU, self.OnTailLength, id=xrc.XRCID("menu_playback_tails") )
        self.frame.Bind( wx.EVT_MENU, self.OnCheckDim, id=xrc.XRCID("menu_playback_dim") )
        self.frame.Bind( wx.EVT_MENU, self.OnCheckZoom, id=xrc.XRCID("menu_settings_zoom") )

        # DEBUG
        if DEBUG:
            self.frame.Bind( wx.EVT_MENU, self.OnDebugPrintIntervals,id=ID_PRINTINTERVALS )
            self.frame.Bind( wx.EVT_MENU, self.OnDebugPrintBufferElement,id=ID_PRINTBUFFEREL )
            self.frame.Bind( wx.EVT_MENU, self.OnDebugPrintFrame,id=ID_PRINTFRAME )

        # stats menu
        #self.frame.Bind( wx.EVT_MENU, self.OnStats, id=xrc.XRCID("menu_stats_vel") )
        #self.frame.Bind( wx.EVT_MENU, self.OnStats, id=xrc.XRCID("menu_stats_orn") )
        #self.frame.Bind( wx.EVT_MENU, self.OnStats, id=xrc.XRCID("menu_stats_space") )
        #self.frame.Bind( wx.EVT_MENU, self.OnStats, id=xrc.XRCID("menu_stats_pos") )
        #self.frame.Bind( wx.EVT_MENU, self.OnSettingsStats, id=xrc.XRCID("menu_stats_settings") )
        #self.frame.Bind( wx.EVT_MENU, self.OnCheckBatchStats, id=xrc.XRCID("menu_stats_batch") )
        # miscellaneous events
        self.frame.Bind( wx.EVT_CLOSE, self.OnQuit )
        self.frame.Bind( wx.EVT_SCROLL, self.OnSlider, self.slider )
        self.frame.Bind( wx.EVT_SIZE, self.OnResize )
        self.frame.Bind( wx.EVT_MAXIMIZE, self.OnResize ) # not called in Gnome?
        # toolbar
        self.frame.Bind(wx.EVT_TOOL,self.ZoomToggle,id=self.zoom_id)
        self.frame.Bind(wx.EVT_TOOL,self.OnPlayButton,id=self.play_id)
        self.frame.Bind(wx.EVT_TOOL,self.OnStopButton,id=self.stop_id)
        self.frame.Bind(wx.EVT_TOOL,self.OnSpeedUpButton,id=self.speedup_id)
        self.frame.Bind(wx.EVT_TOOL,self.OnSlowDownButton,id=self.slowdown_id)
        self.frame.Bind(wx.EVT_TOOL,self.OnRefreshButton,id=self.refresh_id)
        wxvt.setup_validated_integer_callback(self.framenumber_text,
                                              xrc.XRCID('text_framenumber'),
                                              self.OnFrameNumberValidated,
                                              pending_color=params.status_blue)

        # default settings
        if sys.platform == 'win32':
            homedir = os.path.join(os.environ['HOMEDRIVE'], os.environ['HOMEPATH'])
        else:
            homedir = os.environ['HOME']
        self.dir = homedir
        self.ellipse_thickness = params.ellipse_thickness

        # read saved settings, overwriting defaults
        self.ReadUserfile()

        self.OnResize(None)

    def InitState(self):
        # initialize state values
        self.start_frame = 0
        self.tracking = False
        self.batch = None
        self.batch_data = None
        self.bg_window_open = False
        params.batch_executing = False
        params.framesbetweenrefresh = 1
        self.dowritesbfmf = False

    def SetPlayToolTip(self,s):
        self.toolbar.SetToolShortHelp(self.play_id,s)

    def SetStopToolTip(self,s):
        self.toolbar.SetToolShortHelp(self.stop_id,s)

    def SetSpeedUpToolTip(self,s):
        self.toolbar.SetToolShortHelp(self.speedup_id,s)

    def SetSlowDownToolTip(self,s):
        self.toolbar.SetToolShortHelp(self.slowdown_id,s)

    def EnableRefreshBitmap(self,state):
        self.toolbar.EnableTool(self.refresh_id,state)

    def EnablePlayBitmap(self,state):
        self.toolbar.EnableTool(self.play_id,state)

    def EnableStopBitmap(self,state):
        self.toolbar.EnableTool(self.stop_id,state)

    def UpdateToolBar(self,state):
        if state == 'stopped':
            self.EnablePlayBitmap(True)
            self.SetPlayToolTip(self.start_playback_tooltip)
            self.EnableStopBitmap(True)
            self.SetSpeedUpToolTip(self.speedup_playback_tooltip)
            self.SetSlowDownToolTip(self.slowdown_playback_tooltip)
            self.EnableRefreshBitmap(False)
            self.rate_text.SetValue("Play Speed: %.1f fps"%self.play_speed)
        else:
            self.EnableStopBitmap(True)
            self.EnablePlayBitmap(False)
            if self.tracking:
                self.SetStopToolTip(self.stop_tracking_tooltip)
                self.SetSpeedUpToolTip(self.speedup_tracking_tooltip)
                self.SetSlowDownToolTip(self.slowdown_tracking_tooltip)
                self.EnableRefreshBitmap(True)
                if params.do_refresh:
                    self.rate_text.SetValue("Refresh Period: %02d fr"%params.framesbetweenrefresh)
                else:
                    self.rate_text.SetValue("Refresh Rate: Never")
            else:
                self.SetStopToolTip(self.stop_playback_tooltip)
                self.SetSpeedUpToolTip(self.speedup_playback_tooltip)
                self.SetSlowDownToolTip(self.slowdown_playback_tooltip)
                self.EnableRefreshBitmap(False)
                self.rate_text.SetValue("Play Speed: %.1f fps"%self.play_speed)

    def ReadUserfile( self ):
        """Read user settings from file. Set window size, location.
        Called on startup."""
        try:
            if sys.platform[:3] == 'win':
                homedir = os.path.join(os.environ['HOMEDRIVE'], os.environ['HOMEPATH'])
            else:
                homedir = os.environ['HOME']
            self.userfile_name = os.path.join(homedir, '.Ctraxrc' )
            userfile = open( self.userfile_name, "r" )
            # line 1: last version used
            last_version = userfile.readline().rstrip()
            if last_version != __version__:
                wx.MessageBox( "Welcome to Ctrax version %s"%__version__,
                               "Ctrax updated", wx.ICON_INFORMATION )
            # line 2: data directory
            self.dir = userfile.readline().rstrip()
            # line 3: window position and size
            last_pos, last_size = str222tuples( userfile.readline().rstrip() )
            self.frame.SetPosition( last_pos )
            self.frame.SetSize( last_size )
            # line 4: bg window info
            try:
                self.last_bg_pos, self.last_bg_size = \
                                  str222tuples( userfile.readline().rstrip() )
            except (ValueError, IndexError):
                self.last_bg_pos = None
                self.last_bg_size = None
            # line 5: batch window info
            try:
                self.last_batch_pos, self.last_batch_size = \
                                  str222tuples( userfile.readline().rstrip() )
            except (ValueError, IndexError):
                self.last_batch_pos = None
                self.last_batch_size = None
            # line 6: use colorblind pallette?
            use_cb = userfile.readline().rstrip()
            #if use_cb == "True":
            #    self.menu.Check( xrc.XRCID("menu_settings_use_colorblind"), True )
            #    params.use_colorblind_palette = True
            #    params.colors = params.colorblind_palette
            # line 7: ellipse thickness
            use_et = userfile.readline().rstrip()
            if use_et.isdigit():
                self.ellipse_thickness = int(use_et)
            # line 8: tail length
            use_tl = userfile.readline().rstrip()
            if use_tl.isdigit():
                params.tail_length = int(use_tl)
            # line 9: dim image?
            use_di = userfile.readline().rstrip()
            if use_di == "True":
                self.menu.Check( xrc.XRCID("menu_playback_dim"), True )
            # line 10: file-save directory
            self.save_dir = userfile.readline().rstrip()
            draw.const.save_dir = self.save_dir
            # line 11: zoom window info
            try:
                self.last_zoom_pos, self.last_zoom_size = \
                                  str222tuples( userfile.readline().rstrip() )
            except (ValueError, IndexError):
                self.last_zoom_pos = None
                self.last_zoom_size = None
            # done
            userfile.close()
        except:
            # silent on all errors... not really important if this fails
            # however, some attributes do need to be initialized, or
            # else we'll raise errors later (on quit, at latest)
            if not hasattr( self, 'last_bg_pos' ):
                self.last_bg_pos = None
                self.last_bg_size = None
            if not hasattr( self, 'last_batch_pos' ):
                self.last_batch_pos = None
                self.last_batch_size = None
            if not hasattr( self, 'save_dir' ):
                if hasattr( self, 'dir' ):
                    self.save_dir = self.dir
                else:
                    self.save_dir = os.environ['HOME'] # Linux only
            if not hasattr( self, 'last_zoom_pos' ):
                self.last_zoom_pos = None
                self.last_zoom_size = None

    def WriteUserfile( self ):
        """Write current user settings to file. Called on quit."""
        try:
            userfile = open( self.userfile_name, "w" )
        except IOError: pass
        else:
            print >> userfile, __version__ # last app version used
            print >> userfile, self.dir # working directory
            print >> userfile, \
                  self.frame.GetPosition(), \
                  self.frame.GetSize() # window info
            if self.last_bg_pos is not None:
                print >> userfile, \
                      self.last_bg_pos, self.last_bg_size # bg window info
            else: print >> userfile, ""
            if self.last_batch_pos is not None:
                print >> userfile, \
                      self.last_batch_pos, self.last_batch_size # batch window info
            else: print >> userfile, ""
            print >> userfile, False
            print >> userfile, self.ellipse_thickness
            print >> userfile, params.tail_length
            print >> userfile, self.menu.IsChecked( xrc.XRCID("menu_playback_dim") )
            print >> userfile, self.save_dir # file-save directory
            if self.last_zoom_pos is not None:
                print >> userfile, \
                      self.last_zoom_pos, self.last_zoom_size # zoom window info
            else: print >> userfile, ""
            userfile.close()

    def OnFrameNumberValidated( self,evt ):
        new_frame = int(self.framenumber_text.GetValue())
        if new_frame < 0:
            new_frame = 0
        elif new_frame >= self.n_frames:
            new_frame = self.n_frames - 1
        self.start_frame = new_frame
        self.ShowCurrentFrame()

    def OnSettingsBG( self, evt ):
        """Open window for bg threshold settings."""
        # don't create a duplicate window
        if self.bg_window_open:
            self.bg_imgs.frame.Raise()
            return

        # grab previously used threshold, if any, so 'reset' button will work
        if self.ann_file is not None:
            old_thresh = params.n_bg_std_thresh
        else: old_thresh = None

        # show warning if background image calculation is necessary
        isbgmodel = self.CheckForBGModel()
        if isbgmodel == False:
            return

        # set up bg window
        self.bg_imgs.ShowBG( self.frame, self.start_frame, old_thresh )
        self.bg_imgs.frame.Bind( wx.EVT_SIZE, self.OnResizeBG )
        self.bg_imgs.frame.Bind( wx.EVT_MOVE, self.OnResizeBG )
        self.bg_imgs.frame.Bind( wx.EVT_MENU, self.OnQuitBG, id=xrc.XRCID("menu_window_close") )
        self.bg_imgs.frame.Bind( wx.EVT_CLOSE, self.OnQuitBG )
        self.bg_imgs.frame.Bind( wx.EVT_BUTTON, self.OnQuitBG, id=xrc.XRCID("done_button") )

        # set size from memory
        if self.last_bg_pos is not None:
            self.bg_imgs.frame.SetPosition( self.last_bg_pos )
            self.bg_imgs.frame.SetSize( self.last_bg_size )

        # update bg window items with new size
        self.bg_imgs.OnThreshSlider( None )
        self.bg_imgs.OnFrameSlider( None )

        # finally, show bg window
        self.bg_imgs.frame.Show()
        self.bg_window_open = True

        # resize slider and image
        self.bg_imgs.DoSub()
        self.bg_imgs.frame_slider.SetMinSize( wx.Size(
            self.bg_imgs.img_panel.GetRect().GetWidth(),
            self.bg_imgs.frame_slider.GetRect().GetHeight() ) )

    def OnSettingsBGModel( self, evt ):
        """Open window for bg model settings."""

        if params.movie.type == 'sbfmf':
            resp = wx.MessageBox( "Background Model is already set for SBFMF files, and cannot be changed", "Cannot Change Background Model", wx.OK )
            return

        # set up bg window
        if not hasattr( self.bg_imgs, 'modeldlg' ):
            self.bg_imgs.modeldlg = bg.BgSettingsDialog( self.frame, self.bg_imgs )
            self.bg_imgs.modeldlg.frame.Bind( wx.EVT_CLOSE, self.OnQuitBGModel )
            self.bg_imgs.modeldlg.frame.Bind( wx.EVT_BUTTON, self.OnQuitBGModel, id=xrc.XRCID("done_button") )
        else:
            self.bg_imgs.modeldlg.frame.Raise()

    def OnSettingsTracking( self, evt ):
        """Open window for bg model settings."""

        # don't create a duplicate window
        if hasattr(self,'tracking_settings_window_open') and self.tracking_settings_window_open:
            self.tracking_settings.frame.Raise()
            return

        isbgmodel = self.CheckForBGModel()

        if isbgmodel == False:
            return

        # create window
        self.tracking_settings = tracking_settings.TrackingSettings(self.frame,self.bg_imgs,self.start_frame)
        self.tracking_settings.frame.Bind(wx.EVT_CLOSE,self.OnQuitTrackingSettings)
        self.tracking_settings.frame.Bind(wx.EVT_BUTTON,self.OnQuitTrackingSettings, id=xrc.XRCID("done") )

        # finally, show bg window
        self.tracking_settings.frame.Show()
        self.tracking_settings_window_open = True

    def OnChooseOrientations( self, evt ):
        """Open window for bg model settings."""

        # don't create a duplicate window
        if hasattr(self,'choose_orientations_window_open') and self.choose_orientations_window_open:
            self.choose_orientations.frame.Raise()
            return

        # create window
        self.choose_orientations = chooseorientations.ChooseOrientations(self.frame)
        self.choose_orientations.frame.Bind(wx.EVT_CLOSE,self.OnQuitChooseOrientations)
        self.choose_orientations.frame.Bind(wx.EVT_BUTTON,self.OnQuitChooseOrientations, id=xrc.XRCID("ID_CANCEL") )
        self.choose_orientations.frame.Bind(wx.EVT_BUTTON,self.ChooseOrientations, id=xrc.XRCID("ID_OK") )

        # finally, show bg window
        self.choose_orientations.frame.Show()
        self.choose_orientations_window_open = True

    def OnQuitBGModel(self, evt):
        self.bg_imgs.modeldlg.frame.Destroy()
        delattr( self.bg_imgs, 'modeldlg' )

    def OnQuitTrackingSettings(self, evt):

        self.tracking_settings.frame.Destroy()
        delattr( self, 'tracking_settings' )
        self.tracking_settings_window_open = False

    def OnQuitChooseOrientations(self, evt):
        self.choose_orientations.frame.Destroy()
        delattr( self, 'choose_orientations' )
        self.choose_orientations_window_open = False
        #self.RewriteTracks()

    def ChooseOrientations(self, evt):
        self.ann_file = self.choose_orientations.ChooseOrientations(self.ann_file)
        self.OnQuitChooseOrientations(evt)

    def OnResizeBG( self, evt ):
        """BG window was moved or resized. Rescale image and slider,
        and remember new location."""
        if evt is not None: evt.Skip()
        self.bg_imgs.frame.Layout()
        try:
            #self.bg_imgs.DoSub()
            self.bg_imgs.DoDraw()
            self.bg_imgs.frame_slider.SetMinSize( wx.Size(
                self.bg_imgs.img_panel.GetRect().GetWidth(),
                self.bg_imgs.frame_slider.GetRect().GetHeight() ) )
            self.last_bg_size = self.bg_imgs.frame.GetSize()
            self.last_bg_pos = self.bg_imgs.frame.GetPosition()
        except AttributeError: pass # during initialization

    def OnQuitBG( self, evt ):
        """Take data from bg threshold window and close it."""
        self.bg_window_open = False
        self.bg_imgs.hf.frame.Destroy()
        self.bg_imgs.frame.Destroy()
        #wx.MessageBox( "Background threshold set to %d"%self.bg_imgs.show_thresh, "New threshold", wx.ICON_INFORMATION )

    def OnCheckZoom( self, evt ):
        """Open ellipse zoom window."""
        if self.menu.IsChecked( xrc.XRCID("menu_settings_zoom") ):
            # open zoom window
            self.zoom_window = ell.EllipseFrame( self.frame )
            self.zoom_window.frame.Bind( wx.EVT_CLOSE, self.OnCloseZoom )
            self.zoom_window.frame.Bind( wx.EVT_SIZE, self.OnZoomResize )
            self.zoom_window.frame.Bind( wx.EVT_MOVE, self.OnZoomMove )

            # set window position from memory
            if self.last_zoom_pos is not None:
                self.zoom_window.frame.SetPosition( self.last_zoom_pos )
                self.zoom_window.frame.SetSize( self.last_zoom_size )

            if evt is not None: self.ShowCurrentFrame()
        else:
            self.OnCloseZoom( None )

    def OnZoomResize( self, evt ):
        """Zoom window was resized; remember position and redraw."""
        evt.Skip()
        self.zoom_window.Redraw()
        try:
            self.last_zoom_size = self.zoom_window.frame.GetSize()
            self.last_zoom_pos = self.zoom_window.frame.GetPosition()
        except AttributeError: pass # during initialization

    def OnZoomMove( self, evt ):
        """Zoom window moved; remember position."""
        evt.Skip()
        try:
            self.last_zoom_size = self.zoom_window.frame.GetSize()
            self.last_zoom_pos = self.zoom_window.frame.GetPosition()
        except AttributeError: pass # during initialization

    def OnCloseZoom( self, evt ):
        """Close ellipse zoom window."""
        self.menu.Check( xrc.XRCID("menu_settings_zoom"), False )
        self.zoom_window.frame.Destroy()

    def ZoomToggle(self,evt):
        self.zoommode = self.zoommode == False

    def MouseClick(self,evt):

        if not self.zoommode or \
                not self.ann_file.IsAnnData or \
                self.start_frame > self.ann_file.lastframetracked or \
                self.start_frame < self.ann_file.firstframetracked:
            return

        # check to see if zoom window open

        # get the clicked position
        windowheight = self.img_wind_child.GetRect().GetHeight()
        windowwidth = self.img_wind_child.GetRect().GetWidth()
        x = evt.GetX() * self.img_size[1] / windowwidth
        y = self.img_size[0] - evt.GetY() * self.img_size[0] / windowheight
        if (x > self.img_size[1]) or (y > self.img_size[0]):
            return

        # determine closest target
        mind = num.inf
            
        ells = self.ann_file.get_frame(self.start_frame)
        for i,v in ells.iteritems():
            d = (v.center.x-x)**2 + (v.center.y-y)**2
            if d <= mind:
                mini = i
                mind = d
        maxdshowinfo = (num.maximum(params.movie_size[0],params.movie_size[1])/params.MAXDSHOWINFO)**2

        if mind < maxdshowinfo:
            window = None
            if not self.menu.IsChecked( xrc.XRCID("menu_settings_zoom") ):
                # open if not open
                self.menu.Check(xrc.XRCID("menu_settings_zoom"),True)
                self.OnCheckZoom(evt)
                window = 0

            self.ZoomTarget(ells[mini],window)

    def ZoomTarget(self,targ,window):

        # check to see if the target is already drawn in some window
        nwindows = self.zoom_window.n_ell
        for i in range(nwindows):
            if targ.identity == self.zoom_window.ellipse_windows[i].spinner.GetValue():
                # if it is, then do nothing
                return

        # which window will we draw in?
        maxnwindows = self.zoom_window.n_ell_spinner.GetMax()

        if window is None:
            # can we open another window?
            if nwindows < maxnwindows:
                # then open and draw in a new window
                self.zoom_window.AddEllipseWindow(targ.identity)
            else:
                # if we can't, need to choose a window to replace
                if hasattr(self,'firstzoomwindowcreated'):
                    window = self.firstzoomwindowcreated
                    self.firstzoomwindowcreated = (self.firstzoomwindowcreated+1)%maxnwindows
                else:
                    window = 0
                    self.firstzoomwindowcreated = 1%maxnwindows
                # set the identity for the chosen window
                self.zoom_window.ellipse_windows[window].spinner.SetValue(targ.identity)
                self.zoom_window.ellipse_windows[window].redraw()
        else:
            self.zoom_window.ellipse_windows[window].spinner.SetValue(targ.identity)
            self.zoom_window.ellipse_windows[window].redraw()

    #def OnSettingsStats( self, evt ):
    #    """Open statistics settings dialog."""
    #    if not hasattr( self, 'draw_dialog' ):
    #        self.draw_dialog = draw.DrawSettingsDialog( self.frame )
    #        self.draw_dialog.frame.Bind( wx.EVT_CLOSE, self.OnQuitSettings )
    #    else:
    #        self.draw_dialog.frame.Raise()

    #def OnQuitSettings( self, evt ):
    #    """Close stats settings dialog and enable a new one to be opened."""
    #    self.draw_dialog.frame.Destroy()
    #    delattr( self, 'draw_dialog' )

    #def OnCheckColorblind( self, evt ):
    #    """Colorblind-friendly palette was checked or unchecked."""
    #    if self.menu.IsChecked( xrc.XRCID("menu_settings_use_colorblind") ):
    #        params.use_colorblind_palette = True
    #        params.colors = params.colorblind_palette
    #    else:
    #        params.use_colorblind_palette = False
    #        params.colors = params.normal_palette
    #    # rewrite color lists in annotation data already read
    #    if self.ann_data is not None:
    #        for frame in self.ann_data:
    #            if params.use_colorblind_palette:
    #                frame.colors = params.colorblind_palette
    #            else:
    #                frame.colors = params.normal_palette
    #    self.ShowCurrentFrame()

    def OnCheckDim( self, evt ): self.ShowCurrentFrame()

    #def OnEllipseSize( self, evt ):
    #    """Ellipse size dialog."""
    #    try:
    #        dlg = wx.NumberEntryDialog( self.frame,
    #                                    "Enter new ellipse line thickness",
    #                                    "(0-5 pixels)",
    #                                    "Ellipse Thickness",
    #                                    value=self.ellipse_thickness,
    #                                    min=0,
    #                                    max=5 )
    #        if dlg.ShowModal() == wx.ID_OK:
    #            self.ellipse_thickness = dlg.GetValue()
    #        dlg.Destroy()
    #    except AttributeError: # NumberEntryDialog not present yet in 2.6.3.2
    #        import warnings
    #        warnings.filterwarnings( 'ignore', '', DeprecationWarning )
    #        # but for some reason, GetNumberFromUser is already deprecated
    #        new_num = wx.GetNumberFromUser( "Enter new ellipse line thickness",
    #                                        "(0-5 pixels)",
    #                                        "Ellipse Thickness",
    #                                        self.ellipse_thickness,
    #                                        min=0,
    #                                        max=5,
    #                                        parent=self.frame )
    #        warnings.resetwarnings()
    #        if new_num >= 0:
    #            self.ellipse_thickness = new_num
    #    self.ShowCurrentFrame()

    def OnTailLength( self, evt ):
        try:
            dlg = wx.NumberEntryDialog( self.frame,
                                        "Enter new tail length",
                                        "(0-200 frames)",
                                        "Tail Length",
                                        value=params.tail_length,
                                        min=0,
                                        max=200 )
            if dlg.ShowModal() == wx.ID_OK:
                params.tail_length = dlg.GetValue()
            dlg.Destroy()
        except AttributeError: # NumberEntryDialog not present yet in 2.6.3.2
            import warnings
            warnings.filterwarnings( 'ignore', '', DeprecationWarning )
            # but for some reason, GetNumberFromUser is already deprecated
            new_num = wx.GetNumberFromUser( "Enter new tail length",
                                            "(0-200 frames)",
                                            "Tail Length",
                                            params.tail_length,
                                            min=0,
                                            max=200,
                                            parent=self.frame )
            warnings.resetwarnings()
            if new_num >= 0:
                params.tail_length = new_num
        self.ShowCurrentFrame()


    def IsBGModel(self):
        return hasattr(self.bg_imgs,'center')

    def CheckForBGModel(self):

        # already have the background?
        if hasattr(self.bg_imgs,'center'):
            return True

        # in non-interactive mode, automatically do the computation
        if params.interactive == False:
            # do the calculation
            self.bg_imgs.est_bg()
            return True

        # ask what to do
        if params.use_median:
            algtxt = 'Median'
        else:
            algtxt = 'Mean'
        msgtxt = 'Background model has not been calculated.\nCalculate now using the following parameters?\n\nAlgorithm: %s\nNumber of Frames: %d' %(algtxt,params.n_bg_frames)
        if wx.MessageBox( msgtxt, "Calculate?", wx.OK|wx.CANCEL ) == wx.CANCEL:
            # Don't do the computation
            return False

        # set up for running background calculation
        if self.status is not None:
            start_color = self.status.GetBackgroundColour()
            self.status.SetBackgroundColour( params.status_green )
            self.status.SetStatusText( "calculating background", params.status_box )
        wx.BeginBusyCursor()
        wx.Yield()

        # do the calculation
        self.bg_imgs.est_bg(self.frame)

        # return to normal
        if self.status is not None:
            self.status.SetBackgroundColour( start_color )
            self.status.SetStatusText('',params.status_box)
        wx.EndBusyCursor()

        return True

    def CheckForShapeModel(self):

        # already have computed the shape model
        if params.have_computed_shape:
            return True

        haveshape = num.isinf(params.maxshape.area) == False

        # have read in the shape model
        if (params.movie_name == params.annotation_movie_name) and haveshape:
            return True

        # in non-interactive mode, automatically do the computation
        if params.interactive == False:
            # do the calculation
            ell.est_shape( self.bg_imgs )
            return True

        # ask what to do
        if haveshape:
            msgtxt = 'Shape model has not been automatically computed for this movie. Currently:\n\nMin Area = %.2f\nMax Area = %.2f\n\nDo you want to automatically compute now with the following parameters?\n\nNumber of Frames: %d\nNumber of Standard Deviations: %.2f'%(params.minshape.area,params.maxshape.area,params.n_frames_size,params.n_std_thresh)
        else:
            msgtxt = 'Shape is currently unbounded. Do you want to automatically compute now with the following parameters?\n\nNumber of Frames: %d\nNumber of Standard Deviations: %.2f'%(params.n_frames_size,params.n_std_thresh)
        resp = wx.MessageBox( msgtxt, "Calculate?", wx.YES_NO|wx.CANCEL )
        if resp == wx.NO:
            return True
        elif resp == wx.CANCEL:
            return False

        # set up for running shape calculation
        if self.status is not None:
            start_color = self.status.GetBackgroundColour()
            self.status.SetBackgroundColour( params.status_blue )
            self.status.SetStatusText( "calculating shape", params.status_box )
        wx.BeginBusyCursor()
        wx.Yield()

        # do the calculation
        if params.interactive:
            ell.est_shape(self.bg_imgs,self.frame)
        else:
            ell.est_shape(self.bg_imgs)

        # return to normal
        if self.status is not None:
            self.status.SetBackgroundColour( start_color )
            self.status.SetStatusText('',params.status_box)
        wx.EndBusyCursor()

        return True

    def OnDebugPrintIntervals(self,evt):
        print "DEBUG INTERVALS: framestracked = [%d,%d], framesbuffered = [%d,%d], frameswritten = [%d,%d]"%(self.ann_file.firstframetracked,self.ann_file.lastframetracked,self.ann_file.firstframebuffered,self.ann_file.lastframebuffered,self.ann_file.firstframewritten,self.ann_file.lastframewritten)
        sys.stdout.flush()

    def OnDebugPrintBufferElement(self,evt):
        q = wx.TextEntryDialog(self.frame,"Buffer has %d elements, corresponding to frames %d through %d"%(len(self.ann_file.buffer),self.ann_file.firstframebuffered,self.ann_file.lastframebuffered),"Choose buffer element","0")
        q.ShowModal()
        v = q.GetValue()
        try:
            i = int(v)
            if i < 0 or i >= len(self.ann_file.buffer):
                raise NotImplementedError

        except:
            print "DEBUG BUFFER: Must enter integer between 0 and %d"%(len(self.ann_file.buffer)-1)
        else:
            print "DEBUG BUFFER: buffer[%d] (frame %d) = "%(i,i+self.ann_file.firstframebuffered) + str(self.ann_file.buffer[i])
        sys.stdout.flush()

    def OnDebugPrintFrame(self,evt):
        q = wx.TextEntryDialog(self.frame,"Frames tracked: %d through %d ([%d,%d] in buffer)"%(self.ann_file.firstframetracked,self.ann_file.lastframetracked,self.ann_file.firstframebuffered,self.ann_file.lastframebuffered),"Choose frame to print","0")
        q.ShowModal()
        v = q.GetValue()
        try:
            i = int(v)
            if i < self.ann_file.firstframetracked or \
                    i > self.ann_file.lastframetracked:
                raise NotImplementedError
        except:
            print "DEBUG FRAME: Must enter integer between %d and %d"%(self.ann_file.firstframetracked,self.ann_file.lastframetracked)
        else:
            print "DEBUG FRAME: Frame %d = "%i + str(self.ann_file.get_frame(i))
        sys.stdout.flush()



def str222tuples( string ):
    """Converts string into two 2-tuples."""
    vals = string.split()
    for vv in range( len(vals) ):
        vals[vv] = int( vals[vv].strip( '(), ' ) )
    return (vals[0],vals[1]), (vals[2],vals[3])

def str22tuple( string ):
    """Converts string into a 2-tuple."""
    vals = string.split()
    for vv in range( len(vals) ):
        vals[vv] = int( vals[vv].strip( '(), ' ) )
    return (vals[0],vals[1])

