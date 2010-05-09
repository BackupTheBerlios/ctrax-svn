; setup2.nsi
;
; This script is based on example2.nsi.
;
; It will install Ctrax-script.exe into a directory that the user selects,

;--------------------------------

; The name of the installer
Name "Ctrax-0.1.5.6"

; The file to write
OutFile "Ctrax-0.1.5.6-installer.exe"
Icon 'icons\Ctraxicon.ico'

; The default installation directory
InstallDir "$PROGRAMFILES\Ctrax-0.1"

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\Ctrax-0.1.5.6" "Install_Dir"

; Request application privileges for Windows Vista
RequestExecutionLevel admin

;--------------------------------

; Pages

Page components
Page directory
Page instfiles

UninstPage uninstConfirm
UninstPage instfiles

;--------------------------------

; The stuff to install
Section "Ctrax-0.1.5.6 (required)"

  SectionIn RO
  
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR

  ; Whether to compress files
  SetCompress Auto

  ; Whether to overwrite files
  SetOverwrite IfNewer
  
  ; Put file there
  File "dist\*.*"

  SetOutPath $INSTDIR\icons
  SetCompress Auto
  SetOverwrite IfNewer
  File "dist\icons\*.*"

  SetOutPath $INSTDIR\xrc
  SetCompress Auto
  SetOverwrite IfNewer
  File "dist\xrc\*.*"
  
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\Ctrax-0.1.5.6 "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.6" "DisplayName" "Ctrax-0.1.5.6"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.6" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.6" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.6" "NoRepair" 1
  WriteUninstaller "uninstall.exe"
  
SectionEnd

; Optional section (can be disabled by the user)
Section "Start Menu Shortcuts"

  CreateDirectory "$SMPROGRAMS\Ctrax"
  CreateShortCut "$SMPROGRAMS\Ctrax\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  CreateShortCut "$SMPROGRAMS\Ctrax\Ctrax-0.1.lnk" "$INSTDIR\Ctrax-script.exe" "" "$INSTDIR\Ctraxicon.ico" 0
  
SectionEnd

; Optional section (can be disabled by the user)
Section "Desktop Shortcut"

  CreateShortCut "$DESKTOP\Ctrax-0.1.lnk" "$INSTDIR\Ctrax-script.exe" "" "$INSTDIR\icons\Ctraxicon.ico" 0
  
SectionEnd


;--------------------------------

; Uninstaller

Section "Uninstall"
  
  ; Remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.6"
  DeleteRegKey HKLM SOFTWARE\Ctrax-0.1.5.6

  ; Remove files and uninstaller
  Delete $INSTDIR\Ctrax-script.exe
  Delete $INSTDIR\uninstall.exe

  ; Remove shortcuts, if any
  Delete "$SMPROGRAMS\Ctrax\*.*"
  Delete "$DESKTOP\Ctrax-0.1.lnk"

  ; Remove directories used
  RMDir "$SMPROGRAMS\Ctrax"
  RMDir "$INSTDIR"

SectionEnd
