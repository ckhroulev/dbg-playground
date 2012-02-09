from distutils.core import setup, Extension
import numpy

module1 = Extension('dbgenerator', sources=['dbgenerator.c'],
                    include_dirs=[numpy.get_include()+"/numpy"])

setup(name = 'dbgenerator',
        version='1.0',
        description='Drainage basin generating tools',
        ext_modules = [module1])
