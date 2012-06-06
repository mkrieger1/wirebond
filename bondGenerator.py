#!/usr/bin/python
# -*- coding: utf-8 -*-

import random
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
    rings_ = [bond.ring for bond in pair.bonds]
    # rings 4, 5, 6, 7 are individual signals --> extra spacing
    # rings 0, 1, 2, 3 are power/ground nets --> no extra spacing
    if any(ring in [6, 7] for ring in rings_):
        pair.set_min_dist_pboard(470)
        pair.set_min_dist_pboard_wire(200)
    elif any(ring in [4, 5] for ring in rings_):
        pair.set_min_dist_pboard(470)
        pair.set_min_dist_pboard_wire(200)
    # bonds on different rings
    elif not rings_[0] == rings_[1]:
        pair.set_min_dist_pboard(400)
        pair.set_min_dist_pboard_wire(200)
    # bonds on same ring
    else:
        pair.set_min_dist_pboard(90)
        pair.set_min_dist_pboard_wire(90)

    pair.set_min_dist_pchip_wire(90)

#for bond in bonds:
#    f = bond.p1 - Point2D(-1000, 2000)
#    bond.add_force(f.normed()*bond.l)

# start iteration
start = time.time()

NITER = 5000
DAMP = 0.2
fdist = open('dist_violation.txt', 'w')
for i in range(NITER):
    # repulsion of pads on pcb
    for pair in pairs:
        pair.repulsion_pboard(damp=1.0*DAMP)
        pair.repulsion_pboard_wire(damp=0.75*DAMP)
        pair.repulsion_pchip_wire(damp=1.0*DAMP)

    #mark_for_update(pairs)

    s = '%4i %6.1f %6.1f' % (i,
      max(pair._min_dist_pboard-abs(pair._dist_pboard()) for pair in pairs),
      max([pair._min_dist_pboard_wire-abs(pair._dist_pboard_wire(0)[0])
           for pair in pairs if abs(pair._dist_pboard_wire(0)[1]-0.5) < 0.5] +
          [pair._min_dist_pboard_wire-abs(pair._dist_pboard_wire(1)[0])
           for pair in pairs if abs(pair._dist_pboard_wire(1)[1]-0.5) < 0.5])
      )
    print s; print >> fdist, s

    # add custom forces
    for bond in bonds:
        if bond.padnumber in range(0, 49):
            bond.add_force(Point2D(-0.5, 0))
        elif bond.padnumber in range(49, 98):
            bond.add_force(Point2D(0, -0.2))
        elif bond.padnumber in range(98, 147):
            bond.add_force(Point2D(0.2, 0))
        elif bond.padnumber in range(147, 196):
            bond.add_force(Point2D(0, 0.2))

    # experimental
    for bond in bonds:
        ftot = sum(bond.forces)
        frandmax = abs(ftot)*0.2
        frandx = 2*frandmax*(random.random()-0.5)
        frandy = 2*frandmax*(random.random()-0.5)
        frand = Point2D(frandx, frandy)
        bond.add_force(frand)

    if not i % 10:
        with open('animation/bondspstest_%04i.ps' % i, 'w') as f:
            bonds_output_postscript(bonds, f, withforces=True,
                                    scaleforce=1.0/DAMP)

    for bond in bonds:
        bond.apply_force()
        
stop = time.time()
print 'took %.2f seconds.' % (stop-start)

for ring in rings:
    print list(sorted(int(round(bond.length))
                      for bond in bonds
                      if bond.ring == ring))


with open('bondspstest.ps', 'w') as f:
    bonds_output_postscript(bonds, f)

with open('spadic10_revA.pas', 'w') as f:
    bonds_output_altium(bonds, 'SPADIC10_revA', f)

