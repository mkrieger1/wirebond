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


# get all possible pairs of bonds: Npairs = (Nbonds**2-Nbonds)/2
pairs = create_all_bondpairs(bonds)
print len(pairs), "pairs"

# select all pairs where the two endpoints (p2) are potentially closer together
# than min_distance
pairs_p2 = [p for p in pairs if p.in_range_p2()]
print len(pairs_p2), "in range"


for i in range(100):
    process_all_bonds(bonds, pairs, pairs_p2)



with open('bondspstest.ps', 'w') as f:
    bonds_output_postscript(bonds, f)

#with open('spadic10.pas', 'w') as f:
#    bonds_output_altium(bonds, 'spadic10', f)

