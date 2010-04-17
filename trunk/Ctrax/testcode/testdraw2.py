import matplotlib
matplotlib.use( 'WXAgg' )

import wx
from wx import xrc
#import motmot.wxglvideo.simple_overlay as wxvideo
import wx.glcanvas as wxvideo
import pkg_resources # part of setuptools
import numpy as nx
import movies
import imagesk
RSRC_FILE = pkg_resources.resource_filename( __name__, "test.xrc" )

class TestDrawApp(wx.App):
    def OnInit( self ):
        
        rsrc = xrc.XmlResource(RSRC_FILE)
        self.frame = rsrc.LoadFrame(None,"FRAME")
        self.frame.Show()
        self.img_panel = xrc.XRCCTRL(self.frame,"PANEL")
        box = wx.BoxSizer( wx.VERTICAL )
        self.img_panel.SetSizer( box )
        self.img_wind = wxvideo.GLCanvas( self.img_panel, -1 )
        box.Add( self.img_wind, 1, wx.EXPAND )
        self.img_panel.SetAutoLayout( True )
        self.img_panel.Layout()

        wx.EVT_LEFT_DOWN(self.img_wind,self.MouseClick)
        #self.img_wind.Bind(wx.EVT_LEFT_DOWN,self.MouseClick)

        #self.filename = '/home/kristin/FLIES/data/walking_arena/movie20071009_155327.sbfmf'
        #self.movie = movies.Movie(self.filename,True)
        #imd,stamp = self.movie.get_frame( 1 )
        #print 'imd = ' + str(imd)
        #print 'minv = ' + str(nx.min(imd))
        #im8 = imagesk.double2mono8(imd)
        #print 'im8 = ' + str(im8)
        #print 'minv = ' + str(nx.min(im8))
        #print 'im8.shape = ' + str(im8.shape)

        #im8 = nx.zeros((1024,1024),dtype=nx.uint8)

        #self.img_wind.update_image_and_drawings('camera',im8,format='MONO8')

        return True

    def MouseClick(self,evt):
        
        print 'mouse clicked'

def main():
    app = TestDrawApp( 0 )
    app.MainLoop()

if __name__ == '__main__':
    main()
