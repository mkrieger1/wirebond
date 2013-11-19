#!/usr/bin/env python

from distutils.core import setup, Extension
from wirebond import __version__

geom2d_module = Extension('geom2d', sources=['src/geom2d.c',
                                             'src/geom2dmodule.c'])

setup(name='wirebond',
      version=__version__,
      description='Tools for creating wire bond fanouts',
      author='Michael Krieger',
      packages=['wirebond'],
      ext_modules = [geom2d_module]
     )

