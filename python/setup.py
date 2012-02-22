from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

extension = Extension("dbg",
                      sources=["dbg.pyx",
                               "../src/upslope_area.cc",
                               "../src/accumulated_flow.cc",
                               "../src/initialize_mask.cc",
                               "../src/DEM.cc",
                               ],
                      include_dirs=[numpy.get_include(), '../src', '/opt/local/include'],
                      library_dirs=['/opt/local/lib'],
                      libraries=['gsl', 'gslcblas', 'gomp'],
                      extra_compile_args=["-O3", "-ffast-math", "-fopenmp"],
                      language="c++")

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [extension]
      )
