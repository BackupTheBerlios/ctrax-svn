#!/usr/bin/env python

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy
import os, glob
import sys

# include numpy directories for extensions
numpyincludedirs = numpy.get_include()

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
Ctrax_package_data = [ f[6:] for f in glob.glob(os.path.join('Ctrax','xrc','*.xrc'))]+\
                     [ f[6:] for f in glob.glob(os.path.join('Ctrax','icons','*.ico'))]+\
                     [ f[6:] for f in glob.glob(os.path.join('Ctrax','xrc','*.bmp'))]

long_description = """
Ctrax: The Caltech Multiple Fly Tracker

(c) 2007-2010 The Caltech Ethomics Project
http://www.dickinson.caltech.edu/ctrax
branson@caltech.edu

Ctrax is an open-source, freely available, machine vision program for
estimating the positions and orientations of many walking flies,
maintaining their individual identities over long periods of time. It
was designed to allow high-throughput, quantitative analysis of
behavior in freely moving flies. Our primary goal in this project is
to provide quantitative behavior analysis tools to the neuroethology
community, thus we've endeavored to make the system adaptable to other
lab's setups. We have assessed the quality of the tracking results for
our setup, and found that it can maintain fly identities indefinitely
with minimal supervision, and on average for 1.5 fly-hours
automatically.
"""

classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Developers',
    'Intended Audience :: End Users/Desktop',
    'Intended Audience :: Information Technology',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License (GPL)',
    'Natural Language :: English',
    'Operating System :: Microsoft :: Windows',
    'Operating System :: POSIX :: Linux',
    'Programming Language :: C',
    'Programming Language :: C++',
    'Programming Language :: Python :: 2.5',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Artificial Intelligence',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Image Recognition',
    'Topic :: Scientific/Engineering :: Information Analysis',
    'Topic :: Scientific/Engineering :: Medical Science Apps.',
    ]

requires=['cython',
          'motmot.imops',
          'motmot.ufmf',
          'motmot.wxglvideo',
          'motmot.wxvalidatedtext',
          'motmot.wxvideo',
          'numpy',
          'PIL',
          'pygarrayimage',
          'pyglet',
          'scipy',
          'wx',
          ]

setup( 
    name="Ctrax",
    version=this_version,
    author="Caltech Ethomics Project",
    author_email="bransonk@janelia.hhmi.org",
    maintainer="Kristin Branson",
    maintainer_email="bransonk@janelia.hhmi.org",
    url="http://www.dickinson.caltech.edu/Ctrax",
    description="Ctrax: The Caltech Multiple Fly Tracker",
    long_description=long_description,
    download_url="https://developer.berlios.de/projects/ctrax/",
    classifiers=classifiers,
    platforms=['Windows','Linux',],
    packages=['Ctrax'],
    requires=requires,
    provides=['Ctrax',],
    obsoletes=['mtrax',],
    scripts=['Ctrax/Ctrax-script.py','Ctrax/Ctrax.py',],
    cmdclass = {'build_ext': build_ext},
    package_dir={'Ctrax': 'Ctrax'},
    package_data = {'Ctrax':Ctrax_package_data},
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
