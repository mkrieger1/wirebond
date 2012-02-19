#!/usr/bin/python
# -*- coding: utf-8 -*-

from bondGeneratorLib import *
from bondOutputPostscript import *
from bondOutputAltium import *

# read input file
(rings, angles, d_min, bonds) = \
  read_chip_pad_definitions('input/spadic10_pads.txt')

# testing bond manipulating algorithms
#pitch = 95
#(bonds, d_av, d_min) = \
#  iterate_bonds(rings, angles, pitch, bonds, 100, shiftInterposer=False)
#bonds = clean_up(bonds, rings, d_av)

#for bond in bonds:
#    f = bond.p1 - Point2D(-1000, 2000)
#    bond.add_force(f.normed()*bond.l)
#

min_distance = 250

for i in range(100):
    print i
    pairs_done = dict((b, []) for b in bonds)
    for b in bonds:
        neighbors = [n for n in bonds
                     if b.in_range(n, min_distance) and n is not b]
        for n in neighbors:
            if not b in pairs_done[n]:
                b.repulsion(n, min_distance)
                pairs_done[b].append(n)

    for bond in bonds:
        bond.apply_force()


with open('bondspstest.ps', 'w') as f:
    bonds_output_postscript(bonds, f)

#with open('spadic10.pas', 'w') as f:
#    bonds_output_altium(bonds, 'spadic10', f)

