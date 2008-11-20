#!/usr/bin/env python
# Ctraxpostinstall.py
# KMB 11/06/2008

# Ctraxpostinstall.py
import sys
import os.path
import shutil

print "***Running Post-Install Script***"
print "Creating shortcut to Ctrax.exe..."
Ctraxpath = os.path.join(sys.prefix,'Scripts')
description = 'Ctrax: The Caltech Multiple Fly Tracking Program'
iconpath = 'Ctraxicon.ico'
for path in sys.path:
    i = path.find(os.path.join('Ctrax','shortcut'))
    if i >= 0:
        iconpath = os.path.join(path[:i],'Ctrax','Ctraxicon.ico')
create_shortcut("Ctrax.exe",description,"Ctrax.lnk",
                '',Ctraxpath,iconpath)
