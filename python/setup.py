from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

extension = Extension("basins",
                      sources=["basins.pyx", "../src/basins.cc",
                               "../src/streamline.cc",
                               "../src/dem/DEM.cc"],
                      include_dirs=[numpy.get_include(), '../src', '../src/dem'],
                      library_dirs=['/opt/local/lib'],
                      libraries=['gsl', 'gslcblas'],
                      language="c++")

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [extension]
      )
