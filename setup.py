#!/usr/bin/env python

from setuptools import setup, Extension
from wirebond import __version__

geom2d_ext = Extension('wirebond.geom2d',
                       sources=['src/geom2d.c', 'src/geom2d_ext.c'])

wirebond_ext = Extension('wirebond.wirebond_ext',
                         sources=['src/geom2d.c', 'src/wirebond_ext.c'])

setup(name='wirebond',
      version=__version__,
      description='Tools for creating wire bond fanouts',
      author='Michael Krieger',
      packages=['wirebond'],
      ext_modules = [geom2d_ext, wirebond_ext]
     )

