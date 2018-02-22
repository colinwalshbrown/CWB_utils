from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

#ext_modules = [cythonize("cython_util_fns")], ["cython_util_fns.pyx"], include_dirs=[numpy.get_include()])]

setup(
      name = 'Sequencing analysis utility fns',
      ext_modules = cythonize("cython_util_fns.pyx")
      )
