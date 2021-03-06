#!/usr/bin/env python

import numpy
#import numpy.numarray as nn
import os, glob
import sys

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import py2exe

# include directories for hungarian
numpyincludedirs = numpy.get_include()

this_version = '0.1.5.7'

# add all of the .xrc and .bmp files
Ctrax_data_files = [('xrc',glob.glob(os.path.join('xrc','*.xrc'))),
                    ('xrc',glob.glob(os.path.join('xrc','*.bmp'))),
                    ('icons',glob.glob(os.path.join('icons','*.ico')))]
Ctrax_package_data = ['icons/*.ico','xrc/*.xrc','xrc/*.bmp']
print 'Ctrax_package_data: ',Ctrax_package_data
print 'Ctrax_data_files: ',Ctrax_data_files

includes = ['scipy.cluster.vq','scipy.io.matlab.streams',
            'matplotlib.backends.backend_tkagg']

import matplotlib

setup( windows=[{"script": 'Ctrax-script.py',
                 "icon_resources": [(1,"icons/Ctraxicon.ico")]}],
       name="Ctraxexe",
       version=this_version,
       description="The Caltech Multiple Fly Tracker",
       author="Caltech ethomics project",
       author_email="bransonk@janelia.hhmi.org",
       url="http://ctrax.berlios.de",
       cmdclass = {'build_ext': build_ext},
       data_files = Ctrax_data_files,
       package_data = {'Ctrax':Ctrax_package_data},
       ext_modules=[Extension('hungarian',['hungarian.cpp',
                                           'asp.cpp'],
                              include_dirs=[numpyincludedirs,]),
                    Extension('houghcircles_C',
                              ['houghcircles_C.c'],
                              include_dirs=[numpyincludedirs,]),
                    Extension('kcluster2d',
                              ['kcluster2d_cython.pyx'],
                              include_dirs=[numpyincludedirs,]),
                    ],
       options={"py2exe":{"includes":includes}},
       data_files=matplotlib.get_py2exe_datafiles(),
       )
