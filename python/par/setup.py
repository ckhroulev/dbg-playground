from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

extension = Extension("foobar",
                      sources=["foobar.pyx", "foobar_c.c"],
                      include_dirs=[numpy.get_include()],
                      libraries=['gomp'],
                      extra_compile_args=["-O3", "-ffast-math", "-fopenmp"])

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [extension]
      )
