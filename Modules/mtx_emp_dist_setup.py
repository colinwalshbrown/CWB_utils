from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("mtx_emp_dist", ["mtx_emp_dist.pyx"], include_dirs=[numpy.get_include()])]

setup(
      name = 'Distributions for TFBS matrices',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules
      )
