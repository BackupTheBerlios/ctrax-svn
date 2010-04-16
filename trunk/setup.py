#!/usr/bin/env python

from setuptools import setup, Extension
from Cython.Distutils import build_ext
#from distutils.sysconfig import get_python_inc
#from setuptools.dist import Distribution
import numpy
#import numpy.numarray as nn
import os, glob
import sys
#from distutils import version

#from distutils.core import setup
#import py2exe

# include directories for hungarian
numpyincludedirs = numpy.get_include()
#numarrayincludedirs = nn.get_numarray_include_dirs()
#includedirs = numarrayincludedirs+[numpyincludedirs,]

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
       cmdclass = {'build_ext': build_ext},
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
                              include_dirs=[numpyincludedirs,]),
                    Extension('kcluster2d',
                              ['kcluster2d/kcluster2d_cython.pyx'],
                              include_dirs=[numpyincludedirs,]),
                    ]
       )
