from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

import os

prefix = ""
try:
    prefix = os.environ['GSL_PREFIX']
except:
    print "Environment variable GSL_PREFIX not set. Trying known locations..."
    prefixes = ["/usr/", "/usr/local/", "/opt/local/", "/sw/"]

    for path in prefixes:
        print "Checking '%s'..." % path
        try:
            os.stat(path + "include/gsl/gsl_odeiv.h")
            prefix = path
            print "Found GSL in '%s'" % prefix
            break
        except:
            pass

if prefix == "":
    print "Could not find GSL. Stopping..."
    import sys
    sys.exit(1)

extension = Extension("dbg",
                      sources=["python/dbg.pyx",
                               "src/upslope_area.cc",
                               "src/accumulated_flow.cc",
                               "src/initialize_mask.cc",
                               "src/DEM.cc",
                               "src/pnpoly.cc"
                               ],
                      include_dirs=[numpy.get_include(), 'src', prefix + '/include'],
                      library_dirs=[prefix + "/lib"],
                      libraries=['gsl', 'gslcblas', 'gomp'],
                      extra_compile_args=["-O3", "-ffast-math", "-fopenmp", "-Wall"],
                      language="c++")

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [extension]
      )
