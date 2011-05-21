# corresponds to version 0.6.5 of fview and 0.5.8 of flytrax

import shutil
import stat
import os.path
import glob

# directories
import motmot
import pygarrayimage
motmotpath = motmot.__path__
pygarrayimagepath = pygarrayimage.__path__[0]
Ctraxpath = os.path.join('..','..')

# things to copy here with the same base name, use wild cards
tocopy = {}

# Ctrax stuff
tocopy['kcluster2d'] = [os.path.join(Ctraxpath,'kcluster2d','kcluster2d_cython.pyx'),]
tocopy['hungarian'] = glob.glob(os.path.join(Ctraxpath,'hungarian','*.cpp'))+\
    glob.glob(os.path.join(Ctraxpath,'hungarian','*.h'))
tocopy['houghcircles'] = glob.glob(os.path.join(Ctraxpath,'houghcircles','*.c'))
tocopy['Ctrax'] = glob.glob(os.path.join(Ctraxpath,'Ctrax','*.py'))

# motmot stuff
motmot_pkgs = ['imops', 'ufmf', 'wxglvideo', 'wxvalidatedtext', 'wxvideo']
motmot_pkg_files = [['imops.pyd'], ['ufmf.py'], ['wxglvideo.py', 'simple_overlay.py'],
                    ['wxvalidatedtext.py'], ['wxvideo.py']]
for pkg_name, pkg_files in zip( motmot_pkgs, motmot_pkg_files ):
    # find matching motmot path
    for mpath in motmotpath:
        if pkg_name in mpath:
            # append path name for each file needing copy
            pkg_paths = []
            for pfile in pkg_files:
                pkg_paths.append( os.path.join( mpath, pkg_name, pfile ) )
            tocopy[pkg_name] = pkg_paths
    if pkg_name not in tocopy.keys(): print "couldn't find %s"%pkg_name
                
##tocopy['imops'] = [os.path.join(motmotpath,'imops','imops.pyd'),]
##tocopy['ufmf'] = [os.path.join(motmotpath,'ufmf','ufmf.py'),]
##tocopy['wxglvideo'] = [os.path.join(motmotpath,'wxglvideo','wxglvideo.py'),
##                       os.path.join(motmotpath,'wxglvideo','simple_overlay.py'),]
##tocopy['wxvalidatedtext'] = [os.path.join(motmotpath,'wxvalidatedtext','wxvalidatedtext.py'),]
##tocopy['wxvideo'] = [os.path.join(motmotpath,'wxvideo','wxvideo.py'),]

# pygarrayimage
tocopy['pygarrayimage'] = [os.path.join(pygarrayimagepath,'arrayimage.py'),]

# setup script
tocopy['setup_py2exe'] = ['setup_py2exe.py','setup.nsi',]

for (name,srcs) in tocopy.iteritems():
    print '\n%s : '%name
    for src in srcs:
        dst = os.path.join('..',os.path.basename(src))
        print 'cp %s %s'%(src,dst)
        shutil.copy(src,dst)
        os.chmod(dst, stat.S_IRWXU)

def copydir(name,srcdirname,dstdirname,pattern):
    print '\n%s:'%name
    if not os.path.isdir(dstdirname):
        os.mkdir(dstdirname)
        print 'mkdir %s'%dstdirname
    for src in glob.glob(os.path.join(srcdirname,pattern)):
        dst = os.path.join(dstdirname,os.path.basename(src))
        print 'cp %s %s'%(src,dst)
        shutil.copy(src,dst)
        os.chmod(dst, stat.S_IRWXU)
    
# icons directory
srcdirname = os.path.join(Ctraxpath,'Ctrax','icons')
dstdirname = os.path.join('..','icons')
pattern = '*.ico'
name = 'icons'
copydir(name,srcdirname,dstdirname,pattern)

# xrc directory
srcdirname = os.path.join(Ctraxpath,'Ctrax','xrc')
dstdirname = os.path.join('..','xrc')
name = 'xrc'
copydir(name,srcdirname,dstdirname,'*.xrc')
copydir(name,srcdirname,dstdirname,'*.bmp')
