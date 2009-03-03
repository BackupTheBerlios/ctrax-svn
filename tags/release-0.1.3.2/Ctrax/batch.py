# batch.py
# KMB 11/06/2208

import os
import wx
from wx import xrc
from params import params

import pkg_resources # part of setuptools
RSRC_FILE = pkg_resources.resource_filename( __name__, "batch.xrc" )

class BatchWindow:
    def __init__( self, parent, directory, file_list=None ):
        self.file_list = []
        if file_list is not None: self.file_list.append( file_list )
        self.dir = directory
        self.ShowWindow( parent )

    def ShowWindow( self, parent ):
        rsrc = xrc.XmlResource( RSRC_FILE )
        self.frame = rsrc.LoadFrame( parent, "frame_Ctrax_batch" )

        # event bindings
        self.frame.Bind( wx.EVT_BUTTON, self.OnButtonAdd, id=xrc.XRCID("button_add") )
        self.frame.Bind( wx.EVT_BUTTON, self.OnButtonRemove, id=xrc.XRCID("button_remove") )
        self.frame.Bind( wx.EVT_BUTTON, self.OnButtonClose, id=xrc.XRCID("button_close") )
        self.frame.Bind( wx.EVT_CLOSE, self.OnButtonClose )

        # button handles
        self.add_button = xrc.XRCCTRL( self.frame, "button_add" )
        self.remove_button = xrc.XRCCTRL( self.frame, "button_remove" )
        self.execute_button = xrc.XRCCTRL( self.frame, "button_execute" )
        
        # textbox handle
        self.list_box = xrc.XRCCTRL( self.frame, "text_list" )
        self.list_box.Set( self.file_list )

        self.arena_choice = xrc.XRCCTRL(self.frame,"arena_choice")
        self.shape_choice = xrc.XRCCTRL(self.frame,"shape_choice")
        self.bg_model_choice = xrc.XRCCTRL(self.frame,"bg_model_choice")
        self.frame.Bind(wx.EVT_CHOICE,self.OnArenaChoice,id=xrc.XRCID("arena_choice"))
        self.frame.Bind(wx.EVT_CHOICE,self.OnShapeChoice,id=xrc.XRCID("shape_choice"))
        self.frame.Bind(wx.EVT_CHOICE,self.OnBgModelChoice,id=xrc.XRCID("bg_model_choice"))        
        if params.batch_autodetect_arena:
            self.arena_choice.SetSelection(0)
        else:
            self.arena_choice.SetSelection(1)
        if params.batch_autodetect_shape:
            self.shape_choice.SetSelection(0)
        else:
            self.bg_model_choice.SetSelection(1)
        if params.batch_autodetect_bg_model:
            self.bg_model_choice.SetSelection(0)
        else:
            self.bg_model_choice.SetSelection(1)

        self.frame.Show()
        self.is_showing = True
        self.executing = False

    def OnArenaChoice(self,evt):
        if evt is None:
            return
        v = self.arena_choice.GetSelection()
        params.batch_autodetect_arena = v == 0

    def OnShapeChoice(self,evt):
        if evt is None:
            return
        v = self.shape_choice.GetSelection()
        params.batch_autodetect_shape = v == 0
        
    def OnBgModelChoice(self,evt):
        if evt is None:
            return
        v = self.bg_model_choice.GetSelection()
        params.batch_autodetect_bg_model = v == 0

    def OnButtonAdd( self, evt ):
        dlg = wx.FileDialog( self.frame, "Select movie", self.dir, "", "FMF Files|*.fmf|StaticBackgroundFMF Files|*.sbfmf", wx.OPEN )
        
        if dlg.ShowModal() == wx.ID_OK:
            self.dir = dlg.GetDirectory()
            newfile = os.path.join( self.dir, dlg.GetFilename() )

            # check for duplicates
            add_flag = True
            for filename in self.file_list:
                if filename == newfile:
                    wx.MessageBox( "File has already been added,\nnot duplicating", "Duplicate", wx.ICON_WARNING )
                    add_flag = False
                    break
                
            if add_flag:
                self.file_list.append( newfile )
                self.list_box.Set( self.file_list )

        dlg.Destroy()

    def OnButtonRemove( self, evt ):
        for ii in reversed( range( len(self.file_list) ) ):
            if self.list_box.IsSelected( ii ):
                # don't remove currently executing job
                if not self.executing or ii != 0:
                    self.file_list.pop( ii )
        self.list_box.Set( self.file_list )

    def OnButtonClose( self, evt ):
        self.frame.Destroy()
        self.is_showing = False
