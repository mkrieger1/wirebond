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
        pair.set_min_dist_pboard(300)
    elif rings[0] == 3: # == rings[1]
        pair.set_min_dist_pboard(300)
    else:
        pair.set_min_dist_pboard(90)
    pair.set_min_dist_pboard_wire(150)
    pair.set_min_dist_pchip_wire(90)

#for bond in bonds:
#    f = bond.p1 - Point2D(-1000, 2000)
#    bond.add_force(f.normed()*bond.l)

# start iteration
start = time.time()

NITER = 100
DAMP = 0.3
for i in range(NITER):
    # repulsion of pads on pcb
    for pair in pairs:
        pair.repulsion_pboard(damp=DAMP)
        pair.repulsion_pboard_wire(damp=DAMP)
        pair.repulsion_pchip_wire(damp=DAMP)

    #mark_for_update(pairs)

    print '%.1f' % max(pair._min_dist_pboard-abs(pair._dist_pboard()) for pair in pairs)

    if not i % 10:
        with open('animation/bondspstest_%04i.ps' % i, 'w') as f:
            bonds_output_postscript(bonds, f, withforces=True,
                                    scaleforce=1.0/DAMP)

    for bond in bonds:
        bond.apply_force()
        
stop = time.time()
print 'took %.2f seconds.' % (stop-start)


with open('bondspstest.ps', 'w') as f:
    bonds_output_postscript(bonds, f)

