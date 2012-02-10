from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

extension = Extension("basins_py",
                      sources=["basins_py.pyx", "basins.cc", "streamline.cc"],
                      include_dirs=[numpy.get_include()],
                      language="c++")

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [extension]
      )
