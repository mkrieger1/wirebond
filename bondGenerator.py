#!/usr/bin/python
# -*- coding: utf-8 -*-

import time

from bondGeneratorLib import *
from bondInputFile import *
from bondOutputPostscript import *
from bondOutputAltium import *

# read input file and create list of neighboring bond pairs
(bonds, rings) = read_bond_definition('input/spadic10_pins_by_number.txt')
pairs = neighbor_pairs(bonds, rings)
#pairs = all_pairs(bonds)

# set minimum distance for bond position on pcb for each pair
for pair in pairs:
    rings = [bond.ring for bond in pair.bonds]
    # rings 4, 5, 6, 7 are individual signals --> extra spacing
    # rings 0, 1, 2, 3 are power/ground nets --> no extra spacing
    if any(ring in [6, 7] for ring in rings):
        pair.set_min_dist_pboard(600)
    elif any(ring in [4, 5] for ring in rings):
        pair.set_min_dist_pboard(400)
    # bonds on different rings
    elif not rings[0] == rings[1]:
        pair.set_min_dist_pboard(600)
    # bonds on same ring
    else:
        pair.set_min_dist_pboard(90)

    pair.set_min_dist_pboard_wire(300)
    pair.set_min_dist_pchip_wire(90)

#for bond in bonds:
#    f = bond.p1 - Point2D(-1000, 2000)
#    bond.add_force(f.normed()*bond.l)

# start iteration
start = time.time()

NITER = 100
DAMP = 0.5
for i in range(NITER):
    # repulsion of pads on pcb
    for pair in pairs:
        pair.repulsion_pboard(damp=0.5*DAMP)
        pair.repulsion_pboard_wire(damp=DAMP)
        pair.repulsion_pchip_wire(damp=DAMP)

    #mark_for_update(pairs)

    print '%.1f' % max(pair._min_dist_pboard-abs(pair._dist_pboard()) for pair in pairs)

    # add custom forces
    for bond in bonds:
        if bond.padnumber in range(0, 49):
            bond.add_force(Point2D(-5, 0))
        #elif bond.padnumber in range(49, 98):
        #    bond.add_force(Point2D(0, -10))
        #elif bond.padnumber in range(98, 147):
        #    bond.add_force(Point2D(10, 0))
        #elif bond.padnumber in range(147, 196):
        #    bond.add_force(Point2D(0, 10))

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

with open('spadic10_revA.pas', 'w') as f:
    bonds_output_altium(bonds, 'SPADIC10_revA', f)

