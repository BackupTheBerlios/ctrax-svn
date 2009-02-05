;NSIS Script For Ctrax-0.1

!define py2exeOutputDirectory 'dist'

;Background Colors
BGGradient 0000FF 000000 FFFFFF

;Title Of Your Application
Name "Ctrax-0.1"

;Do A CRC Check
CRCCheck On

;Output File Name
OutFile "Ctrax-0.1.exe"
Icon 'Ctraxicon.ico'

;The Default Installation Directory
InstallDir "$PROGRAMFILES\Ctrax-0.1"

;The text to prompt the user to enter a directory
DirText "Please select the folder below"

Section "Install"
  ;Install Files
  SetOutPath $INSTDIR
  SetCompress Auto
  SetOverwrite IfNewer
  File '${py2exeOutputDirectory}\*.*'

  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1" "DisplayName" "Ctrax-0.1 (remove only)"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1" "UninstallString" "$INSTDIR\Uninst.exe"
WriteUninstaller "Uninst.exe"
SectionEnd

Section "Shortcuts"
  ;Add Shortcuts
  CreateShortCut "$SMPROGRAMS\Ctrax-0.1.lnk" "$INSTDIR\Ctrax.exe" "" "$INSTDIR\Ctraxicon.ico" 0
  CreateShortCut "$DESKTOP\Ctrax-0.1.lnk" "$INSTDIR\Ctrax.exe" "" "$INSTDIR\Ctraxicon.ico" 0
SectionEnd

UninstallText "This will uninstall Ctrax-0.1 from your system"

Section Uninstall
  ;Delete Files
  Delete "$INSTDIR\Ctrax-0.1.exe"
  Delete "$DESKTOP\Ctrax-0.1.lnk"

  ;Delete Start Menu Shortcuts
  Delete "$SMPROGRAMS\Ctrax-0.1.lnk"

  ;Delete Uninstaller And Unistall Registry Entries
  Delete "$INSTDIR\Uninst.exe"
  DeleteRegKey HKEY_LOCAL_MACHINE "SOFTWARE\Ctrax-0.1"
  DeleteRegKey HKEY_LOCAL_MACHINE "SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1"
  RMDir "$INSTDIR"
SectionEnd