import wx
from wx import xrc
import motmot.wxvideo.wxvideo as wxvideo
#import motmot.wxglvideo.simple_overlay as wxvideo
import imagesk
import numpy as num
import pkg_resources # part of setuptools
import scipy.ndimage.morphology as morph

RSRC_FILE = pkg_resources.resource_filename( __name__, "fixbg.xrc" )

USEGL = False

class FixBG:

    def __init__(self,parent,bgmodel):
        self.parent = parent
        self.bgcenter = bgmodel.center.copy()
        self.bgdev = bgmodel.dev.copy()
        self.bgcenter0 = self.bgcenter.copy()
        self.bgdev0 = self.bgdev.copy()
        self.bgmodel = bgmodel
        rsrc = xrc.XmlResource(RSRC_FILE)
        self.frame = rsrc.LoadFrame(parent,"fixbg_frame")

        self.InitControlHandles()
        self.InitializeValues()
        self.BindCallbacks()
        self.OnResize()
        self.ShowImage()

    def CheckSave(self):

        if not self.issaved:
            
            msgtxt = 'Some changes to background model may not have been saved. Save now?'
            retval = wx.MessageBox( msgtxt, "Save Changes?", wx.YES_NO|wx.CANCEL )
            if retval == wx.YES:
                self.OnSave()
            return retval != wx.CANCEL

        return True

    def RefreshCallback(self,evt):
        self.ShowImage()

    def control(self,ctrlname):
        return xrc.XRCCTRL(self.frame,ctrlname)

    def InitControlHandles(self):

        self.undo_button = self.control('undo_button')
        self.close_button = self.control('close_button')
        self.cancel_button = self.control('cancel_button')

        self.save_button = self.control('save_button')
        self.quit_button = self.control('quit_button')

        self.displaymode_choice = self.control('displaymode_choice')
        self.showpolygons_checkbox = self.control('showpolygons_checkbox')

        self.undo_button.Enable(False)
        self.close_button.Enable(False)

        self.img_panel = self.control('img_panel')
        box = wx.BoxSizer( wx.VERTICAL )
        self.img_panel.SetSizer( box )
        self.img_wind = wxvideo.DynamicImageCanvas( self.img_panel, -1 )
        self.img_wind.set_resize(True)
        box.Add( self.img_wind, 1, wx.EXPAND )
        self.img_panel.SetAutoLayout( True )
        self.img_panel.Layout()

    def InitializeValues(self):

        wx.BeginBusyCursor()
        wx.Yield()

        imblank = num.zeros(self.bgcenter.shape,dtype=num.uint8)
        self.img_wind.update_image_and_drawings('fixbg',imblank,
                                                format='MONO8')
        self.img_wind_child = self.img_wind.get_child_canvas('fixbg')
        self.polygons = []
        self.currpolygon = num.zeros((0,2))
        (self.nr,self.nc) = self.bgcenter.shape

        # all pixel locations in the image,
        # for checking which pixels are in the polygon
        (self.Y,self.X) = num.mgrid[0:self.nr,0:self.nc]

        # initialize to show center
        self.displaymode_choice.SetSelection(0) 

        # initialize to show polygons
        self.showpolygons_checkbox.SetValue(True)

        self.issaved = True

        wx.EndBusyCursor()

    def BindCallbacks(self):

        # mouse click
        self.img_wind_child.Bind(wx.EVT_LEFT_DOWN,self.MouseClick)
        self.img_wind_child.Bind(wx.EVT_LEFT_DCLICK,self.MouseDoubleClick)

        # buttons
        self.frame.Bind( wx.EVT_SIZE, self.OnResize )
        self.frame.Bind( wx.EVT_CHOICE, self.RefreshCallback, self.displaymode_choice)
        self.frame.Bind( wx.EVT_CHECKBOX, self.RefreshCallback, self.showpolygons_checkbox)
        self.frame.Bind( wx.EVT_BUTTON, self.OnClose, self.close_button )
        self.frame.Bind( wx.EVT_BUTTON, self.OnCancel, self.cancel_button )
        self.frame.Bind( wx.EVT_BUTTON, self.OnUndo, self.undo_button )
        self.frame.Bind( wx.EVT_BUTTON, self.OnSave, self.save_button )

    def OnKeyPress(self,evt):

        # only when in the middle of adding a new polygon
        if self.currpolygon.shape[0] == 0:
            return

        kc = evt.GetKeyCode()
        if kc == wx.WXK_ESCAPE:
            self.OnCancel()
        elif kc == wx.WXK_RETURN:
            self.OnClose()

    def AddPolygon(self):

        wx.BeginBusyCursor()
        wx.Yield()
        self.currpolygon = num.r_[self.currpolygon,num.reshape(self.currpolygon[0,:],(1,2))]
        self.polygons.append(self.currpolygon)
        isin = self.fix_hole(self.bgcenter,self.bgcenter,self.currpolygon)
        isin = self.fix_hole(self.bgdev,self.bgdev,self.currpolygon,isin=isin)
        self.issaved = False
        wx.EndBusyCursor()

    def OnClose(self,evt=None):

        # need at least three points to close
        if self.currpolygon.shape[0] > 2:
            self.AddPolygon()
        else:
            wx.MessageBox( "Not enough points selected to close this polygon",
                           "Select More Points", wx.ICON_ERROR )

        self.currpolygon = num.zeros((0,2))
        self.InAddPolygonMode()
        self.ShowImage()

    def OnSave(self,evt=None):

        self.bgmodel.SetCenter(self.bgcenter)
        self.bgmodel.SetDev(self.bgdev)
        self.bgmodel.SetFixBgPolygons(self.polygons)
        
        self.issaved = True

    def OnCancel(self,evt=None):

        self.currpolygon = num.zeros((0,2))
        self.InAddPolygonMode()
        self.ShowImage()

    def OnUndo(self,evt):

        if self.currpolygon.shape[0] > 0:
            self.currpolygon = num.zeros((0,2))
        elif len(self.polygons) == 0:
            return
        else:
            lastpolygon = self.polygons.pop()
            self.RemovePolygon(lastpolygon)

        self.InAddPolygonMode()
        self.ShowImage()

    def ShowImage(self):

        polygoncolor = (255,0,0)
        currpolygoncolor = (255,255,0)

        if self.displaymode_choice.GetSelection() == 0:
            self.img_shown = imagesk.double2mono8(self.bgcenter,donormalize=True)
        else:
            self.img_shown = imagesk.double2mono8(self.bgdev,donormalize=True)

        if self.currpolygon.shape[0] >= 1:
            points = [[self.currpolygon[0,0],self.currpolygon[0,1]]]
            pointcolors = [(currpolygoncolor[0]/255.,currpolygoncolor[1]/255.,
                            currpolygoncolor[2]/255.)]
        else:
            points = []
            pointcolors = []
        
        lines = []
        for j in range(self.currpolygon.shape[0]-1):
            lines.append([self.currpolygon[j,0],
                          self.currpolygon[j,1],
                          self.currpolygon[j+1,0],
                          self.currpolygon[j+1,1],
                          currpolygoncolor])

        if self.showpolygons_checkbox.GetValue():
            for i in range(len(self.polygons)):
                for j in range(self.polygons[i].shape[0]-1):
                    lines.append([self.polygons[i][j,0],
                                  self.polygons[i][j,1],
                                  self.polygons[i][j+1,0],
                                  self.polygons[i][j+1,1],
                                  polygoncolor])

        (lines,linecolors) = imagesk.separate_linesegs_colors(lines)

        self.img_wind.update_image_and_drawings('fixbg',self.img_shown,
                                                format="MONO8",
                                                linesegs=lines,
                                                lineseg_colors=linecolors,
                                                points=points,
                                                point_colors=pointcolors)

        self.img_wind.Refresh(eraseBackground=False)
        
    def MouseClick(self,evt):

        resize = self.img_wind.get_resize()
        x = evt.GetX() / resize
        y = self.nr - evt.GetY() / resize
        if (x > self.nc) or (y > self.nr):
            return
        self.currpolygon = num.r_[self.currpolygon,num.array([x,y],ndmin=2)]
        self.InAddPointMode()
        self.ShowImage()

    def MouseDoubleClick(self,evt):

        # not in the middle of adding points, then treat as a single click
        if self.currpolygon.shape[0] == 0:
            self.MouseClick(evt)
            return

        resize = self.img_wind.get_resize()
        x = evt.GetX() / resize
        y = self.nr - evt.GetY() / resize
        if (x > self.nc) or (y > self.nr):
            return

        # make sure we're close to the first pixel
        d = num.sqrt((x - self.currpolygon[0,0])**2 + (y - self.currpolygon[0,1])**2)
        if d > 10:
            return

        self.OnClose()

    def InAddPointMode(self):
        self.cancel_button.Enable(True)
        self.close_button.Enable(True)
        self.undo_button.Enable(False)

    def InAddPolygonMode(self):
        self.cancel_button.Enable(False)
        self.close_button.Enable(False)
        self.undo_button.Enable(len(self.polygons) > 0)

    def RemovePolygon(self,poly):
        wx.BeginBusyCursor()
        wx.Yield()
        isin = self.undo_fix_hole(self.bgcenter0,self.bgcenter,poly)
        isin = self.undo_fix_hole(self.bgdev0,self.bgdev,poly,isin=isin)
        self.issaved = False
        wx.EndBusyCursor()

    def fix_hole(self,im,out,polygon,isin=None):

        s = num.ones((3,3),bool)        

        # get points on the inside of the polygon
        if isin is None:
            isin = point_inside_polygon(self.X,self.Y,polygon)
        (y_isin,x_isin) = num.nonzero(isin)

        # get points on the outside boundary of the polygon
        isboundary = num.logical_and(morph.binary_dilation(isin,s),
                                     ~isin)

        (y_isboundary,x_isboundary) = num.nonzero(isboundary)

        # for each point inside the polygon, get an average from
        # all the boundary points
        for i in range(len(y_isin)):
            x = x_isin[i]
            y = y_isin[i]
            d = num.sqrt((y_isboundary-y)**2+(x_isboundary-x)**2)
            w = num.exp(-d)
            out[y,x] = num.sum(im[isboundary] * w) / num.sum(w)

        return isin

    def undo_fix_hole(self,im0,im,polygon,isin=None):
        if isin is None:
            isin = point_inside_polygon(self.X,self.Y,polygon)
        im[isin] = im0[isin]    
        return isin

    def OnResize(self,evt=None):

        if evt is not None: evt.Skip()
        self.frame.Layout()
        try:

            #new_size = wx.Size( self.img_wind.GetRect().GetWidth(),
            #                    self.displaymode_choice.GetRect().GetHeight() )
            #self.displaymode_choice.SetMinSize(new_size)
            #self.displaymode_choice.SetSize(new_size)

            self.ShowImage()

        except AttributeError: pass # during initialization

def point_inside_polygon(x,y,poly):
    # adapted from http://www.ariel.com.au/a/python-point-int-poly.html
    # by Patrick Jordan

    n = len(poly)
    inside = num.zeros(x.shape,dtype=bool)
    xinters = num.zeros(x.shape)

    p1x,p1y = poly[0]
    idx2 = num.zeros(x.shape,dtype=bool)
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        idx = num.logical_and(y > min(p1y,p2y),
                              num.logical_and(y <= max(p1y,p2y),
                                              x <= max(p1x,p2x)))
        if p1y != p2y:
            xinters[idx] = (y[idx]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
        if p1x == p2x:
            inside[idx] = ~inside[idx]
        else:
            idx2[:] = False
            idx2[idx] = x[idx] <= xinters[idx]
            inside[idx2] = ~inside[idx2]
        p1x,p1y = p2x,p2y

    return inside
