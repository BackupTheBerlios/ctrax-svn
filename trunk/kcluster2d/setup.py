from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

numpyincludedirs = numpy.get_include()

ext_modules = [#Extension("convolve1", ["convolve1.pyx"],
               #          include_dirs=[numpyincludedirs,]),
               #Extension("convolve2", ["convolve2.pyx"],
               #          include_dirs=[numpyincludedirs,]),
               Extension("kcluster", ["kcluster.pyx"],
                         include_dirs=[numpyincludedirs,]),
               Extension("loop", ["loop.pyx"],
                         include_dirs=[numpyincludedirs,]),
               ]

setup(
  name = 'test cython extensions',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
