#!/usr/bin/env python

try:
    from setuptools import setup, Extension
except ImportError:
    print 'Installing setuptools ...'
    import ez_setup
    ez_setup.use_setuptools()
    from setuptools import setup, Extension
    print 'Done installing setuptools.'
from distutils.sysconfig import get_python_inc
from setuptools.dist import Distribution
import numpy
import numpy.numarray as nn
import os, glob
import sys
from distutils import version

#from distutils.core import setup
#import py2exe

# include directories for hungarian
numpyincludedirs = numpy.get_include()
#numarrayincludedirs = nn.get_numarray_include_dirs()
#includedirs = numarrayincludedirs+[numpyincludedirs,]

kws = {}
if int(os.getenv( 'DIABLE_INSTALL_REQUIRES','1' )):
    print "Setuptools install_requires disabled"
else:
    # commented wx, numpy, scipy, PIL out because eggs are not 
    # available for them. From now on, DISABLE_INSTALL_REQUIRES
    # is by default 1
    print "Checking requirements using Setuptools install_requires"
    install_requires=[#'wxPython>=2.8',
                      #'numpy>=1.0.3',
                      #'scipy',
                      #'PIL>=1.1.6',
                      'motmot.wxvideo>=0.5.2.dev',
                      'motmot.wxglvideo>=0.6.1.dev',
                      'motmot.wxvalidatedtext',
                      'motmot.imops',
                      'pygarrayimage',
                      'pyglet>=1.0',
                      ]
    kws['install_requires'] = install_requires

# read version number from version file
path = os.path.abspath( os.curdir )
Ctrax_path = os.path.join( path, 'Ctrax' )
ver_filename = os.path.join( Ctrax_path, 'version.py' )
ver_file = open( ver_filename, "r" )
for line in ver_file: # parse through file version.py
    if line.find( '__version__' ) >= 0:
        line_sp = line.split() # split by whitespace
        version_str = line_sp[2] # third item
        this_version = version_str[1:-1] # strip quotes
ver_file.close()

# add all of the .xrc and .bmp files
Ctrax_package_data = [ f[6:] for f in glob.glob(os.path.join('Ctrax','*.xrc'))]+\
                     [ f[6:] for f in glob.glob(os.path.join('Ctrax','*.bmp'))]+\
                     [ 'Ctraxicon.ico']
eager_resources = [ f for f in glob.glob(os.path.join('Ctrax','*.xrc')) ] + \
    [ f for f in glob.glob(os.path.join('Ctrax','*.bmp'))] + \
    [ f for f in glob.glob('Ctraxicon.ico)')]
print 'Ctrax_package_data: ',Ctrax_package_data

setup( name="Ctrax",
       version=this_version,
       description="Multiple fly tracker",
       author="Caltech ethomics project",
       author_email="branson@caltech.edu",
       url="http://www.dickinson.caltech.edu/Ctrax",
       packages=['Ctrax'],
       entry_points = {'console_scripts': ['Ctrax=Ctrax:main']},
       package_dir={'Ctrax': 'Ctrax'},
       #py_modules=['Ctrax.colormapk','Ctrax.imagesk',
       #            'Ctrax.houghcircles','Ctrax.setarena'],
       package_data = {'Ctrax':Ctrax_package_data},
       eager_resources=eager_resources,
       dependency_links = [
           "http://alldunn.com/wxPython/stuff/",
           #"http://vision.caltech.edu/~kristin/pythoneggs_for_Ctrax/",
           ],
       ext_modules=[Extension('hungarian',['hungarian/hungarian.cpp',
                                           'hungarian/asp.cpp'],
                              include_dirs=[numpyincludedirs,]),
                    Extension('houghcircles_C',
                              ['houghcircles/houghcircles_C.c'],
                              include_dirs=[numpyincludedirs,])
                    ],
       **kws)
