#!/usr/bin/env python

import cProfile, pstats
import os
import shutil
import subprocess

# rebuild
try:
    shutil.rmtree('build')
except OSError:
    pass
subprocess.call(['python', 'setup.py', 'build'])

# profile
os.chdir('examples')

pr = cProfile.Profile()
pr.enable()
execfile('spadic10_revA.py')
pr.disable()

stats = pstats.Stats(pr)
stats.strip_dirs() # use before calc_callees()
stats.calc_callees()

def print_calls(func, caller=None, level=0):
    callees = stats.all_callees[func]
    if not caller:
        data = stats.stats[func][:-1]
    else:
        data = stats.all_callees[caller][func]
    print_calls_line(func, data, level)
    for f in callees:
        print_calls(f, func, level+1)

def print_calls_line(func, data, level):
    nc, cc, tt, ct = data
    print (' '*level*7 + "%6.3f %6.3f" % (ct, tt) + 
           ' '*(8-level)*7 + "%8i %s" %(nc, pstats.func_std_string(func)))

# get spadic10_revA.run() call
width, result = stats.get_print_list(['spadic.*run'])
func = result[0]
print_calls(func)

