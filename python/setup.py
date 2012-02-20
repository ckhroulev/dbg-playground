from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

extension = Extension("basins",
                      sources=["basins.pyx",
                               "../src/basins.cc",
                               "../src/streamline.cc",
                               "../src/DEM.cc",
                               "../src/init_mask.cc"],
                      include_dirs=[numpy.get_include(), '../src'],
                      library_dirs=['/opt/local/lib'],
                      libraries=['gsl', 'gslcblas', 'gomp'],
                      extra_compile_args=["-O3", "-ffast-math", "-fopenmp"],
                      language="c++")

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [extension]
      )
