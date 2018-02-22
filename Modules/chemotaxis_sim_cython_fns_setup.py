from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

include_gsl_dir = "/opt/local/include/"
lib_gsl_dir = "/opt/local/lib/"

ext_modules = [Extension("chemotaxis_sim_cython_fns",
                          ["chemotaxis_sim_cython_fns.pyx"],
                          include_dirs=[numpy.get_include(),include_gsl_dir],
                          library_dirs=[lib_gsl_dir],
                          libraries=["gsl"])]

setup(
      name = 'chemotaxis_sim_cython_fns',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules
      )
