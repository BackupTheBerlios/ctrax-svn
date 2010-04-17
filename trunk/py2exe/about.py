import wx

# forward compatibility for AboutDialogInfo and AboutBox, added in wx 2.7.1.1
# (currently 2.6.3.2 is packaged for Ubuntu Feisty)
class AboutDialogInfo:
    def __init__( self ):
	self.name=self.description=self.version=self.copyright=""
    def SetName( self, name ): self.name = name
    def SetDescription( self, desc ): self.description = desc
    def SetVersion( self, version ): self.version = version
    def SetCopyright( self, cprt ): self.copyright = cprt
def AboutBox( info ):
    wx.MessageBox( u"%s version %s\n\n%s\n\ncopyright \u00A9%s"%(info.name, info.version, info.description, info.copyright), "About %s"%info.name, wx.ICON_INFORMATION )

