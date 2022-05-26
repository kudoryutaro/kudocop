from distutils.core import setup
from Cython.Build import cythonize

setup(name="Molcop",
      ext_modules=cythonize("**/*.pyx"))
