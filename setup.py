#!/usr/bin/env python

from distutils.core import setup
from wirebond import __version__

setup(name='wirebond',
      version=__version__,
      description='Tools for creating wire bond fanouts',
      author='Michael Krieger',
      packages=['wirebond'],
     )

