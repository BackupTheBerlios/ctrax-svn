import wx
from wx import xrc
import motmot.wxglvideo.simple_overlay as wxvideo
import pkg_resources # part of setuptools
import numpy as nx
RSRC_FILE = pkg_resources.resource_filename( __name__, "test.xrc" )

class TestDrawApp(wx.App):
    def OnInit( self ):
        
        rsrc = xrc.XmlResource(RSRC_FILE)
        self.frame = rsrc.LoadFrame(None,"FRAME")
        self.frame.Show()
        self.img_panel = xrc.XRCCTRL(self.frame,"PANEL")
        box = wx.BoxSizer( wx.VERTICAL )
        self.img_panel.SetSizer( box )
        self.img_wind = wxvideo.DynamicImageCanvas( self.img_panel, -1 )
        box.Add( self.img_wind, 1, wx.EXPAND )
        self.img_panel.SetAutoLayout( True )
        self.img_panel.Layout()
        print 'here!'
        im8 = nx.zeros((1024,1024),dtype=nx.uint8)

        self.img_wind.update_image_and_drawings('camera',im8,format='MONO8')

        return True

def main():
    app = TestDrawApp( 0 )
    app.MainLoop()

if __name__ == '__main__':
    main()
