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

# get all possible pairs of bonds: Npairs = (Nbonds**2-Nbonds)/2
pairs = []
for i in range(len(bonds)):
    for j in range(i+1, len(bonds)):
        pairs.append(BondPair(bonds[i], bonds[j]))
print len(pairs), "pairs"

# select all pairs where the two endpoints (p2) are potentially closer together
# than min_distance
pairs_in_range = [p for p in pairs if p.in_range_p2(min_distance)]
print len(pairs_in_range), "in range"

for i in range(100):
#---------------------------------------------------------------------
# begin iteration
#---------------------------------------------------------------------
    #-----------------------
    # pair subset 1: p2-p2 distance
    #-----------------------

    # sort pairs_in_range by their p2 distance
    pairs_in_range.sort() # only needed for displaying min. distance
    p = pairs_in_range[0]
    print "min. p2 distance:", str(p), p.dist_p2,

    # add repulsive forces if necessary
    for p in pairs_in_range:
        if p.dist_p2 < min_distance:
            p.repulsion_p2(min_distance)

    #-----------------------
    # pair subset 2: p2-wire distance
    #-----------------------

    #-----------------------
    # update all bonds
    #-----------------------
    for p in pairs_in_range:
        if any(len(b.forces) for b in p.pair):
            p.needs_update = True
    for b in bonds:
        b.apply_force()

    #-----------------------
    # update all pairs
    #-----------------------
    updated = 0
    for p in pairs:
        if p.needs_update:
            p.update()
            updated += 1
    print updated, "pairs updated"

#---------------------------------------------------------------------
# end iteration
#---------------------------------------------------------------------



with open('bondspstest.ps', 'w') as f:
    bonds_output_postscript(bonds, f)

#with open('spadic10.pas', 'w') as f:
#    bonds_output_altium(bonds, 'spadic10', f)

