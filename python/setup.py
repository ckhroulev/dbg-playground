from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

extension = Extension("foo",
                      sources=["foo.pyx", "cfoo.c"],
                      include_dirs=[numpy.get_include()])

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [extension]
      )
