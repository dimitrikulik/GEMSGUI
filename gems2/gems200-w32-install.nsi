; gems200-w32-install.nsi
; requires NullSoft NSIS v.1.98 or later to be installed on your PC!
;
; This script will install GEM-Selektor program package 
; and database into a directory that the user selects

; The name of the installer
Name "GEMS2-PSI"
Caption "GEM-Selektor Package"
BrandingText "(C) 2000-2003, GEMS-PSI Team"
Icon "C:\GEMS2inst\program\img\gems2.ico"
WindowIcon on

; The file to write
; OutFile "gems200-w32-install.exe"
OutFile "gems200rc4-w32-install.exe"

; The default installation directory
InstallDir $PROGRAMFILES\GEMS2
; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM SOFTWARE\GEMS2 "Install_Dir"

;License text
LicenseText "Version to be installed: 2.0.0-PSI $\nLICENSE and DISCLAIMER"
LicenseData "C:\GEMS2inst\program\doc\txt\TermsOfUse.txt"

; The text to prompt the user to enter a directory
ComponentText "The GEM-Selektor v.2-PSI package is about to be installed on your PC. (De)select optional components that you want to install, and hit 'Next' to continue."
; The text to prompt the user to enter a directory
DirText "At best, choose something like D:\MyGEMS2 to make finding your \projects easier; or simply hit 'Install' to proceed (with an upgrade)"

; The stuff to install
Section "GEMS2 (required)"
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR
  ; Put file there
  File "C:\GEMS2inst\readme1.txt"
  File "C:\GEMS2inst\rungems2.bat"
  CreateDirectory "$INSTDIR\program"
  SetOutPath $INSTDIR\program
  Delete "$INSTDIR\program\qt*.dll"
  ; Put file there
    File "C:\GEMS2inst\program\gems2.exe"
    File "C:\GEMS2inst\program\gems2.ico"
    File "C:\GEMS2inst\program\qtmt312.dll"
    File "C:\GEMS2inst\program\cc3250mt.dll"
    File "C:\GEMS2inst\program\imm32.dll"
    File "C:\GEMS2inst\program\readme.txt"
  CreateDirectory "$INSTDIR\program\data"
  SetOutPath $INSTDIR\program\data
    File "C:\GEMS2inst\program\data\units.ini"
    File "C:\GEMS2inst\program\data\vis_cn.ini"
    File "C:\GEMS2inst\program\data\vis_od.ini"
    File "C:\GEMS2inst\program\data\viz_cndoc.txt"
  CreateDirectory "$INSTDIR\program\visor.data"
  SetOutPath "$INSTDIR\program\visor.data"
    File "C:\GEMS2inst\program\visor.data\visor.dat"
    File "C:\GEMS2inst\program\visor.data\visobj.dat"
  CreateDirectory "$INSTDIR\program\img"
  SetOutPath $INSTDIR\program\img
   Delete "$INSTDIR\program\img\*.*"
    File "C:\GEMS2inst\program\img\*.*"
  CreateDirectory "$INSTDIR\program\DB.default"
  SetOutPath "$INSTDIR\program\DB.default"
   Delete "$INSTDIR\program\DB.default\*.*"
    File "C:\GEMS2inst\program\DB.default\*.ndx"
    File "C:\GEMS2inst\program\DB.default\*.pdb"
  CreateDirectory "$INSTDIR\program\doc"
  CreateDirectory "$INSTDIR\program\doc\html"
  CreateDirectory "$INSTDIR\program\doc\pdf"
  CreateDirectory "$INSTDIR\program\doc\txt"
  SetOutPath $INSTDIR\program\doc\txt
   Delete "$INSTDIR\program\doc\txt\*.*" 
    File "C:\GEMS2inst\program\doc\txt\*.txt"
  SetOutPath $INSTDIR\program\doc\html
   Delete "$INSTDIR\program\doc\html\*.*"
    File "C:\GEMS2inst\program\doc\html\*.html"
    File "C:\GEMS2inst\program\doc\html\*.png"
    File "C:\GEMS2inst\program\doc\html\prefDialog.jpg"
  SetOutPath $INSTDIR\program\doc\pdf
   Delete "$INSTDIR\program\doc\pdf\*.*"
    File "C:\GEMS2inst\program\doc\pdf\*.pdf"
  SetOutPath $INSTDIR
  CreateDirectory "$INSTDIR\projects"
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\GEMS2 "Install_Dir" "$INSTDIR"
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\GEMS2" "DisplayName" "GEMS2 (remove only)"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\GEMS2" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteUninstaller "uninstall.exe"
SectionEnd

; optional section
Section "Start Menu Shortcuts"
  CreateDirectory "$SMPROGRAMS\GEM-Selektor v2-PSI"
  CreateShortCut "$SMPROGRAMS\GEM-Selektor v2-PSI\Uninstall GEMS2.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  CreateShortCut "$SMPROGRAMS\GEM-Selektor v2-PSI\Run GEMS2.lnk" "$INSTDIR\program\gems2.exe" '-c -s "$INSTDIR\program" -u "$INSTDIR"' "$INSTDIR\program\gems2.ico" 0
  CreateShortCut "$DESKTOP\Run GEMS2.lnk" "$INSTDIR\program\gems2.exe" '-c -s "$INSTDIR\program" -u "$INSTDIR"' "$INSTDIR\program\gems2.ico" 0
SectionEnd

; optional section
Section "Test Modelling Projects"
  ; Set output path to the installation directory.
;  CreateDirectory "$INSTDIR\projects"
  SetOutPath $INSTDIR\projects
  ; Put files here 
    File "C:\GEMS2inst\projects\*.conf"
  CreateDirectory "$INSTDIR\projects\TryNPTDB"
  SetOutPath $INSTDIR\projects\TryNPTDB
    File "C:\GEMS2inst\projects\TryNPTDB\*.pdb"
    File "C:\GEMS2inst\projects\TryNPTDB\*.ndx"
; More test projects to come here!  
;
  SetOutPath $INSTDIR
SectionEnd

; uninstall stuff
UninstallText "This will uninstall the GEM-Selektor v.2-PSI package. Your modelling project subdirectories will not be removed. Hit 'Uninstall' to continue."

; special uninstall section.
Section "Uninstall"
; remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\GEMS2"
  DeleteRegKey HKLM SOFTWARE\GEMS2
; remove files
  Delete "$INSTDIR\*.*"
  Delete "$INSTDIR\program\data\*.*"
  Delete "$INSTDIR\program\visor.data\*.*"
  Delete "$INSTDIR\program\img\*.*"
  Delete "$INSTDIR\program\DB.default\*.*"
  Delete "$INSTDIR\program\*.*" 
  Delete "$INSTDIR\program\doc\txt\*.*"
  Delete "$INSTDIR\program\doc\html\*.*"
  Delete "$INSTDIR\program\doc\pdf\*.*"
;
  Delete "$INSTDIR\projects\*.*" 
;  Delete "$INSTDIR\projects\TryNPTDB\*.*" 
;
; MUST REMOVE UNINSTALLER, too
  Delete $INSTDIR\uninstall.exe
; remove shortcuts, if any.
  Delete "$SMPROGRAMS\GEM-Selektor\*.*"
  Delete "$DESKTOP\Run GEMS2.lnk"
; remove directories used.
  RMDir "$SMPROGRAMS\GEM-Selektor"
  RMDir "$INSTDIR\program\data"
  RMDir "$INSTDIR\program\visor.data"
  RMDir "$INSTDIR\program\img"
  RMDir "$INSTDIR\program\DB.default"
  RMDir "$INSTDIR\program\doc\html"
  RMDir "$INSTDIR\program\doc\pdf"
  RMDir "$INSTDIR\program\doc\txt"
  RMDir "$INSTDIR\program\doc"
  RMDir "$INSTDIR\program" 
;
  RMDir "$INSTDIR\projects\TryNPTDB" 
;  RMDir "$INSTDIR\projects"
;
;  RMDir "$INSTDIR"
SectionEnd

; eof