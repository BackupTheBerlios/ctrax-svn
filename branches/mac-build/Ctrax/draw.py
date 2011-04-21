# draw.py
# KMB 11/10/08

# matplotlib stripped out
##import matplotlib
##matplotlib.interactive( True )
###matplotlib.use( 'WXAgg' )
##import matplotlib.backends.backend_wxagg
##import matplotlib.figure
##import matplotlib.cm
import numpy as num
import os
import wx
from wx import xrc
import imagesk

import motmot.wxvalidatedtext.wxvalidatedtext as wxvt # part of Motmot

from ellipsesk import draw_ellipses,draw_ellipses_bmp

from params import params

##
##class PlotConstants:
##    def __init__( self ):
##        self.vel_bins = 100
##        self.vel_x_min=self.vel_x_max=self.vel_y_min=self.vel_y_max = 'a'
##        self.orn_bins = 90
##        self.orn_x_min=self.orn_x_max=self.orn_y_min=self.orn_y_max = 'a'
##        self.space_bins = 50
##        self.space_x_min=self.space_x_max=self.space_y_min=self.space_y_max = 'a'
##        self.pos_binsize = 50
##        self.pos_x_min=self.pos_x_max=self.pos_y_min=self.pos_y_max = 'a'
##
##        self.save_dir = None
##const = PlotConstants()

def annotate_image( ellipses=None, old_pts=None, thickness=params.ellipse_thickness ):
    """Draw stuff on image."""

    # draw ellipses
    if ellipses is None:
        linesegs = []
    else:
        linesegs = draw_ellipses( ellipses )

    # draw tails
    if old_pts is not None:

        prevpts = {}
        for frame_pts in old_pts:
            for pt in frame_pts:
                if prevpts.has_key(pt[2]):
                    prevpt = prevpts[pt[2]]
                    color = params.colors[pt[2]%len(params.colors)]
                    linesegs.append([prevpt[0]+1,prevpt[1]+1,
                                     pt[0]+1,pt[1]+1,color])
                prevpts[pt[2]] = pt[0:2]

    return linesegs

def annotate_bmp( img, ellipses=None, old_pts=None, thickness=params.ellipse_thickness, windowsize=None, zoomaxes=None ):
    """Draw stuff on image."""

    # draw ellipses
    (bmp,resize,img_size) = draw_ellipses_bmp( img, ellipses, thickness=thickness,
                                               windowsize=windowsize, zoomaxes=zoomaxes )

    # draw tails
    if old_pts is not None:

        # create list of lines
        linedict = {}
        for frame_pts in old_pts:
            for pt in frame_pts:
                if linedict.has_key(pt[2]):
                    linedict[pt[2]].append([pt[0]+1,pt[1]+1])
                else:
                    linedict[pt[2]] = [[pt[0]+1,pt[1]+1],]
        linelists = []
        linecolors = []
        for i,j in linedict.iteritems():
            linecolors.append(params.colors[i%len(params.colors)])
            linelists.append(j)

        bmp = imagesk.add_annotations(bmp,resize,linelists=linelists,linecolors=linecolors,linewidths=[thickness,])

    return (bmp,resize,img_size)


### class from http://www.scipy.org/Matplotlib_figure_in_a_wx_panel
##class NoRepaintCanvas(matplotlib.backends.backend_wxagg.FigureCanvasWxAgg):
##    """We subclass FigureCanvasWxAgg, overriding the _onPaint method, so that
##    the draw method is only called for the first two paint events. After that,
##    the canvas will only be redrawn when it is resized."""
##    def __init__(self, *args, **kwargs):
##        matplotlib.backends.backend_wxagg.FigureCanvasWxAgg.__init__(self, *args, **kwargs)
##        self._drawn = 0
##
##    def _onPaint(self, evt):
##        """
##        Called when wxPaintEvt is generated
##        """
##        if not self._isRealized:
##            self.realize()
##        if self._drawn < 2:
##            self.draw(repaint = False)
##            self._drawn += 1
##        self.gui_repaint(drawDC=wx.PaintDC(self))
##        
### class from http://www.scipy.org/Matplotlib_figure_in_a_wx_panel
### although menu/right-click/save functionality was added here
##class PlotPanel(wx.Panel):
##    """The PlotPanel has a Figure and a Canvas. OnSize events simply set a 
##    flag, and the actually redrawing of the figure is triggered by an Idle event."""
##    def __init__(self, parent, id = -1, color = None,\
##        dpi = None, style = wx.NO_FULL_REPAINT_ON_RESIZE, **kwargs):
##        wx.Panel.__init__(self, parent, id = id, style = style, **kwargs)
##        self.figure = matplotlib.figure.Figure(None, dpi)
##        self.canvas = NoRepaintCanvas(self, -1, self.figure)
##        self.SetColor(color)
##        self.Bind(wx.EVT_IDLE, self._onIdle)
##        self.Bind(wx.EVT_SIZE, self._onSize)
##        self._resizeflag = True
##        self._SetSize()
##
##        # generate context menu for saving files
##        self.menu = wx.Menu()
##        file_save_item = wx.MenuItem( self.menu, wx.ID_ANY, "Save as..." )
##        self.menu.AppendItem( file_save_item )
##        self.canvas.Bind( wx.EVT_RIGHT_UP, self.OnRightMouseButton )
##        self.canvas.Bind( wx.EVT_MENU, self.OnMenuSave )
##
##    def OnRightMouseButton( self, evt ):
##        """Right mouse button pressed; pop up save menu."""
##        self.canvas.PopupMenu( self.menu )
##
##    def OnMenuSave( self, evt ):
##        """User has chosen to save this figure as an image file.
##        Prompt for filename and save."""
##        # the extension on filename determines file format
##        # create text list of allowed file extensions
##        allowed_ext = ('eps', 'png', 'ps', 'svg') # 'pdf' broken in matplotlib
##        allowed_text = "'%s'"%allowed_ext[0]
##        for ext in allowed_ext[1:]:
##            allowed_text += ", '%s'"%ext
##
##        # show dialog
##        dlg = wx.FileDialog( self.canvas, "Save as image of types %s"%allowed_text, const.save_dir, style=wx.SAVE )
##        if dlg.ShowModal() == wx.ID_OK:
##            # get entered filename and extension
##            this_file = dlg.GetFilename()
##            const.save_dir = dlg.GetDirectory()
##            filename = os.path.join( const.save_dir, this_file )
##            this_name, this_ext = filename.rsplit( '.', 1 )
##            if this_ext in allowed_ext:
##                try:
##                    self.figure.savefig( filename )
##                except IOError:
##                    wx.MessageBox( "Miscellaneous write error",
##                                   "Write error", wx.ICON_ERROR )
##                    # love those miscellaneous errors!
##            else:
##                wx.MessageBox( "Filename extension must be one of\n%s"%allowed_text, "Bad file type", wx.ICON_ERROR )
##
##        dlg.Destroy()
##
##    def SetColor(self, rgbtuple):
##        """Set figure and canvas colours to be the same."""
##        if not rgbtuple:
##            rgbtuple = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE).Get()
##        col = [c/255.0 for c in rgbtuple]
##        self.figure.set_facecolor(col)
##        self.figure.set_edgecolor(col)
##        self.canvas.SetBackgroundColour(wx.Colour(*rgbtuple))
##
##    def _onSize(self, event):
##        self._resizeflag = True
##
##    def _onIdle(self, evt):
##        if self._resizeflag:
##            self._resizeflag = False
##            self._SetSize()
##            self.draw()
##
##    def _SetSize(self, pixels = None):
##        """This method can be called to force the Plot to be a desired size,
##        which defaults to the ClientSize of the panel."""
##        
##        if not pixels:
##            pixels = self.GetClientSize()
##        self.canvas.SetSize(pixels)
##        self.figure.set_size_inches(pixels[0]/self.figure.get_dpi(),
##        pixels[1]/self.figure.get_dpi())
##
##    def draw(self): pass
##
##class DebugPlotPanel( PlotPanel ):
##
##    def imshow(self,im):
##        self.grid = im
##        print self.grid.shape
##        print self.grid.T
##
##    def draw(self):
##        if not hasattr(self, 'subplot'):
##            self.subplot = self.figure.add_subplot(111)
##        if hasattr(self,'grid'):
##            self.subplot.imshow( self.grid, cmap=matplotlib.cm.jet )
##
##class VelPlotPanel( PlotPanel ):
##    def __init__( self, frame, data, title, max_jump ):
##        PlotPanel.__init__( self, frame )
##        self.data = data
##        self.title = title
##
##        self.vels = []
##        for ff in range( len(self.data)-1 ):
##            for fly in self.data[ff].iterkeys():
##                if self.data[ff+1].hasItem(fly):
##                    vel = self.data[ff][fly].Euc_dist( self.data[ff+1][fly] )
##                    if vel < max_jump:
##                        self.vels.append( vel )
##                # calculate velocity
##            # for each fly
##        # for each pair of consecutive frames in data
##        self.vels = num.array( self.vels )
##
##    def draw( self ):
##        if not hasattr(self, 'subplot'):
##            self.subplot = self.figure.add_subplot(111)
##
##        # histogram
##        self.subplot.hist( self.vels, const.vel_bins )
##        
##        #Set some plot attributes
##        self.subplot.set_title(self.title, fontsize = 12)
##        self.subplot.set_xlabel( "fly velocity (pixels/frame)" )
##        self.subplot.set_ylabel( "bin count" )
##        if const.vel_x_min != 'a':
##            self.subplot.axis( xmin=const.vel_x_min )
##        if const.vel_x_max != 'a':
##            self.subplot.axis( xmax=const.vel_x_max )
##        if const.vel_y_min != 'a':
##            self.subplot.axis( ymin=const.vel_y_min )
##        if const.vel_y_max != 'a':
##            self.subplot.axis( ymax=const.vel_y_max )
##
##class OrnPlotPanel( PlotPanel ):
##    def __init__( self, frame, data, title ):
##        PlotPanel.__init__( self, frame )
##        self.data = data
##        self.title = title
##
##        self.orns = []
##        for frame in self.data:
##            for fly in frame.itervalues():
##                self.orns.append( fly.angle )
##                # grab angle
##            # for each ellipse
##        # for each frame in data
##        self.orns = num.array( self.orns )
##        
##    def draw( self ):
##        if not hasattr(self, 'subplot'):
##            self.subplot = self.figure.add_subplot(111)
##
##        # histogram
##        self.subplot.hist( self.orns, const.orn_bins )
##        
##        #Set some plot attributes
##        self.subplot.set_title(self.title, fontsize = 12)
##        self.subplot.set_xlabel( "fly orientation (rad)" )
##        self.subplot.set_ylabel( "bin count" )
##        
##class SpacePlotPanel( PlotPanel ):
##    def __init__( self, frame, data, title ):
##        PlotPanel.__init__( self, frame )
##        self.data = data
##        self.title = title
##        
##        self.dists = []
##        for frame in self.data:
##            for fly1 in frame.iterkeys():
##                for fly2 in frame.iterkeys():
##                    if fly1 < fly2:
##                        self.dists.append( frame[fly1].Euc_dist( frame[fly2] ) )
##                    # calculate distance
##                # for each pair of ellipses
##        # for each frame in data
##        self.dists = num.array( self.dists )
##
##    def draw( self ):
##        if not hasattr(self, 'subplot'):
##            self.subplot = self.figure.add_subplot(111)
##
##        # histogram
##        self.subplot.hist( self.dists, const.space_bins )
##        
##        #Set some plot attributes
##        self.subplot.set_title(self.title, fontsize = 12)
##        self.subplot.set_xlabel( "distance between flies (pixels)" )
##        self.subplot.set_ylabel( "bin count" )
##        
##class PosPlotPanel( PlotPanel ):
##    def __init__( self, frame, data, title, width, height ):
##        PlotPanel.__init__( self, frame )
##        self.data = data
##        self.title = title
##        self.width = width
##        self.height = height
##
##        # fill grid
##        self.grid = num.zeros( (self.width/const.pos_binsize+1, self.height/const.pos_binsize+1) )
##        for frame in self.data:
##            for fly in frame.itervalues():
##                fly_x = fly.center.x
##                fly_y = fly.center.y
##                try: self.grid[int(fly_x/const.pos_binsize),int(fly_y/const.pos_binsize)] += 1
##                except IndexError:
##                    print "error adding", fly_x, fly_y, "to bin", fly_x/const.pos_binsize, fly_y/const.pos_binsize, "(size is", self.grid.shape, "\b)"
##                # increment grid
##            # for each fly
##        # for each frame
##
##        self.xticks = range( 0, self.grid.shape[0], max( 200/const.pos_binsize, 1 ) )
##        self.xticklabels = []
##        for xx in self.xticks:
##            self.xticklabels.append( str(xx*const.pos_binsize) )
##        self.yticks = range( 0, self.grid.shape[1], max( 200/const.pos_binsize, 1 ) )
##        self.yticklabels = []
##        for yy in self.yticks:
##            self.yticklabels.append( str(yy*const.pos_binsize) )
##
##    def draw( self ):
##        if not hasattr(self, 'subplot'):
##            self.subplot = self.figure.add_subplot(111)
##
##        # 2D histogram
##        self.subplot.imshow( self.grid.T, cmap=matplotlib.cm.hot )
##
##        #Set some plot attributes
##        self.subplot.set_title(self.title, fontsize = 12)
##        self.subplot.set_xlabel( "x (pixels)" )
##        self.subplot.set_ylabel( "y (pixels)" )
##        self.subplot.set_xticks( self.xticks )
##        self.subplot.set_yticks( self.yticks )
##        self.subplot.set_xticklabels( self.xticklabels )
##        self.subplot.set_yticklabels( self.yticklabels )
##
##class DrawSettingsDialog:
##    def __init__( self, parent=None ):
##        rsrc = xrc.XmlResource( RSRC_FILE )
##        self.frame = rsrc.LoadFrame( parent, "frame_Ctrax_settings" )
##
##        # input box handles
##        self.vel_bins_box = xrc.XRCCTRL( self.frame, "text_vel_bins" )
##        self.orn_bins_box = xrc.XRCCTRL( self.frame, "text_orn_bins" )
##        self.space_bins_box = xrc.XRCCTRL( self.frame, "text_space_bins" )
##        self.pos_binsize_box = xrc.XRCCTRL( self.frame, "text_pos_binsize" )
##        self.vel_axes_box = xrc.XRCCTRL( self.frame, "text_vel_axes" )
##        self.orn_axes_box = xrc.XRCCTRL( self.frame, "text_orn_axes" )
##        self.space_axes_box = xrc.XRCCTRL( self.frame, "text_space_axes" )
##        self.pos_axes_box = xrc.XRCCTRL( self.frame, "text_pos_axes" )
##
##        # initial values
##        self.vel_bins_box.SetValue( str(const.vel_bins) )
##        self.orn_bins_box.SetValue( str(const.orn_bins) )
##        self.space_bins_box.SetValue( str(const.space_bins) )
##        self.pos_binsize_box.SetValue( str(const.pos_binsize) )
##        self.vel_axes_box.SetValue( "a,a,a,a" )
##        self.orn_axes_box.SetValue( "a,a,a,a" )
##        self.space_axes_box.SetValue( "a,a,a,a" )
##        self.pos_axes_box.SetValue( "a,a,a,a" )
##
##        # bind to wxvalidatedtext
##        wxvt.setup_validated_integer_callback( self.vel_bins_box,
##                                               xrc.XRCID("text_vel_bins"),
##                                               self.OnTextValidated,
##                                               pending_color=params.wxvt_bg )
##        wxvt.setup_validated_integer_callback( self.orn_bins_box,
##                                               xrc.XRCID("text_orn_bins"),
##                                               self.OnTextValidated,
##                                               pending_color=params.wxvt_bg )
##        wxvt.setup_validated_integer_callback( self.space_bins_box,
##                                               xrc.XRCID("text_space_bins"),
##                                               self.OnTextValidated,
##                                               pending_color=params.wxvt_bg )
##        wxvt.setup_validated_integer_callback( self.pos_binsize_box,
##                                               xrc.XRCID("text_pos_binsize"),
##                                               self.OnTextValidated,
##                                               pending_color=params.wxvt_bg )
##        wxvt.Validator( self.vel_axes_box, xrc.XRCID("text_vel_axes"),
##                        self.OnTextValidated, self.ParseAxes,
##                        pending_color=params.wxvt_bg )
##        wxvt.Validator( self.orn_axes_box, xrc.XRCID("text_orn_axes"),
##                        self.OnTextValidated, self.ParseAxes,
##                        pending_color=params.wxvt_bg )
##        wxvt.Validator( self.space_axes_box, xrc.XRCID("text_space_axes"),
##                        self.OnTextValidated, self.ParseAxes,
##                        pending_color=params.wxvt_bg )
##        wxvt.Validator( self.pos_axes_box, xrc.XRCID("text_pos_axes"),
##                        self.OnTextValidated, self.ParseAxes,
##                        pending_color=params.wxvt_bg )
##
##        self.frame.Show()
##
##    def OnTextValidated( self, evt ):
##        """User-entered text has been validated; set global variables
##        correspondingly."""
##        if evt.GetId() == xrc.XRCID( "text_vel_bins" ):
##            const.vel_bins = int(self.vel_bins_box.GetValue())
##        elif evt.GetId() == xrc.XRCID( "text_orn_bins" ):
##            const.orn_bins = int(self.orn_bins_box.GetValue())
##        elif evt.GetId() == xrc.XRCID( "text_space_bins" ):
##            const.space_bins = int(self.space_bins_box.GetValue())
##        elif evt.GetId() == xrc.XRCID( "text_pos_binsize" ):
##            const.pos_binsize = int(self.pos_binsize_box.GetValue())
##        elif evt.GetId() == xrc.XRCID( "text_vel_axes" ):
##            const.vel_x_min = self.new_axes_vals[0]
##            const.vel_x_max = self.new_axes_vals[1]
##            const.vel_y_min = self.new_axes_vals[2]
##            const.vel_y_max = self.new_axes_vals[3]
##            # set axes box values to currently stored values (removes user formatting)
##            self.vel_axes_box.SetValue( "%s,%s,%s,%s"%(str(const.vel_x_min),
##                                                       str(const.vel_x_max),
##                                                       str(const.vel_y_min),
##                                                       str(const.vel_y_max)) )
##        elif evt.GetId() == xrc.XRCID( "text_orn_axes" ):
##            const.orn_x_min = self.new_axes_vals[0]
##            const.orn_x_max = self.new_axes_vals[1]
##            const.orn_y_min = self.new_axes_vals[2]
##            const.orn_y_max = self.new_axes_vals[3]
##            # set axes box values to currently stored values (removes user formatting)
##            self.orn_axes_box.SetValue( "%s,%s,%s,%s"%(str(const.orn_x_min),
##                                                       str(const.orn_x_max),
##                                                       str(const.orn_y_min),
##                                                       str(const.orn_y_max)) )
##        elif evt.GetId() == xrc.XRCID( "text_space_axes" ):
##            const.space_x_min = self.new_axes_vals[0]
##            const.space_x_max = self.new_axes_vals[1]
##            const.space_y_min = self.new_axes_vals[2]
##            const.space_y_max = self.new_axes_vals[3]
##            # set axes box values to currently stored values (removes user formatting)
##            self.space_axes_box.SetValue( "%s,%s,%s,%s"%(str(const.space_x_min),
##                                                       str(const.space_x_max),
##                                                       str(const.space_y_min),
##                                                       str(const.space_y_max)) )
##        elif evt.GetId() == xrc.XRCID( "text_pos_axes" ):
##            const.pos_x_min = self.new_axes_vals[0]
##            const.pos_x_max = self.new_axes_vals[1]
##            const.pos_y_min = self.new_axes_vals[2]
##            const.pos_y_max = self.new_axes_vals[3]
##            # set axes box values to currently stored values (removes user formatting)
##            self.pos_axes_box.SetValue( "%s,%s,%s,%s"%(str(const.pos_x_min),
##                                                       str(const.pos_x_max),
##                                                       str(const.pos_y_min),
##                                                       str(const.pos_y_max)) )
##
##    def ParseAxes( self, string ):
##        """Parse user input for drawings axes limits. Return True or False
##        depending on whether input is valid."""
##        vals = string.split( ',' ) # split on commas
##        if len(vals) != 4: vals = string.split() # try again with whitespace
##        if len(vals) != 4:
##            return False
##
##        # strip parentheses, etc. and make numeric
##        for vv in range( len(vals) ):
##            vals[vv] = vals[vv].strip( '()[], ' )
##            # test for empty
##            if vals[vv] == '': return False
##            # test for a
##            if vals[vv] == 'a': continue
##            # test for int
##            if vals[vv].isdigit() or \
##                   (vals[vv][0] == '-' and vals[vv][1:].isdigit()): # negative int
##                vals[vv] = int(vals[vv])
##                continue
##            # test for float
##            try:
##                vals[vv] = float(vals[vv])
##            except ValueError: pass
##            else: continue
##
##            return False
##
##        # validation successful; save values so we don't have to re-parse later
##        self.new_axes_vals = vals
##        return True
##
