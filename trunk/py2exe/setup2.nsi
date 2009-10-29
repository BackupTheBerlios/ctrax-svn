; setup2.nsi
;
; This script is based on example2.nsi.
;
; It will install Ctrax.exe into a directory that the user selects,

;--------------------------------

; The name of the installer
Name "Ctrax-0.1.5.2"

; The file to write
OutFile "Ctrax-0.1.5.2-installer.exe"
Icon 'Ctraxicon.ico'

; The default installation directory
InstallDir "$PROGRAMFILES\Ctrax-0.1"

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\Ctrax-0.1.5.2" "Install_Dir"

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
Section "Ctrax-0.1.5.2 (required)"

  SectionIn RO
  
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR

  ; Whether to compress files
  SetCompress Auto

  ; Whether to overwrite files
  SetOverwrite IfNewer
  
  ; Put file there
  File "dist\*.*"
  
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\Ctrax-0.1.5.2 "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.2" "DisplayName" "Ctrax-0.1.5.2"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.2" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.2" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.2" "NoRepair" 1
  WriteUninstaller "uninstall.exe"
  
SectionEnd

; Optional section (can be disabled by the user)
Section "Start Menu Shortcuts"

  CreateDirectory "$SMPROGRAMS\Ctrax"
  CreateShortCut "$SMPROGRAMS\Ctrax\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  CreateShortCut "$SMPROGRAMS\Ctrax\Ctrax-0.1.lnk" "$INSTDIR\Ctrax.exe" "" "$INSTDIR\Ctraxicon.ico" 0
  
SectionEnd

; Optional section (can be disabled by the user)
Section "Desktop Shortcut"

  CreateShortCut "$DESKTOP\Ctrax-0.1.lnk" "$INSTDIR\Ctrax.exe" "" "$INSTDIR\Ctraxicon.ico" 0
  
SectionEnd


;--------------------------------

; Uninstaller

Section "Uninstall"
  
  ; Remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Ctrax-0.1.5.2"
  DeleteRegKey HKLM SOFTWARE\Ctrax-0.1.5.2

  ; Remove files and uninstaller
  Delete $INSTDIR\Ctrax.exe
  Delete $INSTDIR\uninstall.exe

  ; Remove shortcuts, if any
  Delete "$SMPROGRAMS\Ctrax\*.*"
  Delete "$DESKTOP\Ctrax-0.1.lnk"

  ; Remove directories used
  RMDir "$SMPROGRAMS\Ctrax"
  RMDir "$INSTDIR"

SectionEnd
