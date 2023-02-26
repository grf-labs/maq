import numpy
import os
import sys

from setuptools import setup, find_packages
from distutils.core import Extension
from Cython.Build import cythonize

if 'darwin' in sys.platform:
  COMPILE_ARGS = ['-std=c++11', '-stdlib=libc++', '-Wall', '-O2', '-pthread']
  LINK_ARGS = ['-stdlib=libc++']
elif 'linux' in sys.platform:
  COMPILE_ARGS = ['-std=c++11', '-lstdc++', '-Wall', '-O2', '-pthread']
  LINK_ARGS = ['-lstdc++', '-pthread']
elif 'win32' in sys.platform:
  COMPILE_ARGS = ['-std=c++11', '-lstdc++', '-Wall', '-O2']
  LINK_ARGS = ['-lstdc++']
else
  raise ImportError('Unsupported OS.')

setup_dir = os.path.abspath(os.path.dirname(__file__))
INCLUDE_DIRS = [os.path.join(setup_dir, '..', 'core', 'src')]
INCLUDE_DIRS.append(os.path.join(setup_dir, '..', 'core', 'third_party'))
INCLUDE_DIRS.append(numpy.get_include())
SOURCES = ['maq' + os.path.sep + 'maqbindings.pyx']

ext = Extension(
  'maq.ext',
  language='c++',
  sources=SOURCES,
  extra_compile_args=COMPILE_ARGS,
  extra_link_args=LINK_ARGS,
  include_dirs=INCLUDE_DIRS,
  define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
)

setup(
  name='maq',
  version='0.1.0',
  packages=find_packages(include=['maq']),
  ext_modules=cythonize(ext, compiler_directives={'language_level': 3}),
  url = 'https://github.com/grf-labs/maq'
)
