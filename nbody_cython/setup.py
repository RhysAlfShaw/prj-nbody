from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np
import os

os.environ['NPY_NO_DEPRECATED_API'] = '1'

ext_modules = [
    Extension("nbody", ["nbody_cython/nbody.pyx"],
              include_dirs=[np.get_include()],
              extra_compile_args=["-O3", "-march=native"],
              language="c++")]

setup(
    name="Simulation",
    ext_modules=cythonize(ext_modules, annotate=True),
)