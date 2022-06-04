from distutils.core import setup
from Cython.Build import cythonize

setup(name="KudoCop",
      ext_modules=cythonize("**/*.pyx"))
