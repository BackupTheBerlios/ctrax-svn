#!/usr/bin/env python

import numpy
#import numpy.numarray as nn
import os, glob
import sys

from distutils.core import setup, Extension
import py2exe

# include directories for hungarian
numpyincludedirs = numpy.get_include()

this_version = '0.1'

# add all of the .xrc and .bmp files
Ctrax_package_data = [ f[6:] for f in glob.glob('*.xrc')]+\
                     [ f[6:] for f in glob.glob('*.bmp')]+\
                     [ 'Ctraxicon.ico']
#eager_resources = [ f for f in glob.glob('*.xrc') ] + \
#    [ f for f in glob.glob('*.bmp')] + \
#    [ f for f in glob.glob('Ctraxicon.ico)')]
print 'Ctrax_package_data: ',Ctrax_package_data

setup( console=[{"script": 'Ctrax.py',
                 "icon_resources": [(1,"Ctraxicon.ico")]
                 }
                ],
       name="Ctraxexe",
       version=this_version,
       description="The Caltech Multiple Fly Tracker",
       author="Caltech ethomics project",
       author_email="branson@caltech.edu",
       url="http://www.dickinson.caltech.edu/Research/Ctrax",
       #entry_points = {'console_scripts': ['Ctrax=Ctrax:main']},
       package_data = {'Ctrax':Ctrax_package_data},
       #eager_resources=eager_resources,
       ext_modules=[Extension('hungarian',['hungarian.cpp',
                                           'asp.cpp'],
                              include_dirs=[numpyincludedirs,]),
                    Extension('houghcircles_C',
                              ['houghcircles_C.c'],
                              include_dirs=[numpyincludedirs,])
                    ],
       )
