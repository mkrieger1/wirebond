#!/usr/bin/env python

import os
import shutil
import subprocess

# rebuild
try:
    shutil.rmtree('build')
except OSError:
    pass
subprocess.call(['python', 'setup.py', 'build'])

# run example program
os.chdir('examples')
execfile('spadic10_revA.py')

