#!/usr/bin/python
# -*- coding: utf-8 -*-

import time

from bondGeneratorLib import *
from bondInputFile import *
from bondOutputPostscript import *
from bondOutputAltium import *

# read input file and create list of neighboring bond pairs
(bonds, rings) = read_bond_definition('input/spadic10_pads.txt')
pairs = neighbor_pairs(bonds, rings)
#pairs = all_pairs(bonds)

# set minimum distance for bond position on pcb for each pair
for pair in pairs:
    rings = [bond.ring for bond in pair.bonds]
    if not rings[0] == rings[1]:
        pair.min_dist_pboard = 300
    elif rings[0] == 3: # == rings[1]
        pair.min_dist_pboard = 300
    else:
        pair.min_dist_pboard = 90

#for bond in bonds:
#    f = bond.p1 - Point2D(-1000, 2000)
#    bond.add_force(f.normed()*bond.l)

# start iteration
start = time.time()

NITER = 1000
for i in range(NITER):
    # repulsion of pads on pcb
    for pair in pairs:
        pair.repulsion_pboard(damp=0.3)
        pair.repulsion_pboard_wire(damp=0.3)

    mark_for_update(pairs)

    print '%.1f' % max(pair.min_dist_pboard-pair.dist_pboard() for pair in pairs)

    for bond in bonds:
        bond.apply_force()

    if not i % 100:
        with open('animation/bondspstest_%04i.ps' % i, 'w') as f:
            bonds_output_postscript(bonds, f)
        
stop = time.time()
print 'took %.2f seconds.' % (stop-start)


with open('bondspstest.ps', 'w') as f:
    bonds_output_postscript(bonds, f)

