#!/usr/bin/python
# -*- coding: utf-8 -*-

import random
import time

from wirebond.core import *
from wirebond.input_file import *
from wirebond.output_postscript import *
from wirebond.output_altium import *


class BondGenerator():
#--------------------------------------------------------------------
# setup
#--------------------------------------------------------------------
    def __init__(self, inputfile):
        # read input file and create list of neighboring bond pairs
        (bonds, groups) = read_bond_definition(inputfile)
        pairs = neighbor_pairs(bonds, groups)
        #pairs = all_pairs(bonds)

        # set minimum distance for bond position on pcb for each pair
        for pair in pairs:
            groups_ = [bond.group for bond in pair.bonds]
            # groups 4, 5, 6, 7 are individual signals --> extra spacing
            # groups 0, 1, 2, 3 are power/ground nets --> no extra spacing
            if all(group == 6 for group in groups_):
                pair.set_min_dist_pboard(500)
                pair.set_min_dist_pboard_wire(220)
            elif all(group == 4 for group in groups_):
                pair.set_min_dist_pboard(550)
                pair.set_min_dist_pboard_wire(220)
            elif any(group in [4, 5, 6, 7] for group in groups_):
                pair.set_min_dist_pboard(400)
                pair.set_min_dist_pboard_wire(220)
            # bonds in different groups
            elif not groups_[0] == groups_[1]:
                pair.set_min_dist_pboard(400)
                pair.set_min_dist_pboard_wire(220)
            # bonds in same group
            else:
                pair.set_min_dist_pboard(90)
                pair.set_min_dist_pboard_wire(90)

            pair.set_min_dist_pchip_wire(78)

        self.bonds = bonds
        self.pairs = pairs
        self.groups = groups
        self.it = 0
        self.damp = 0.5

#--------------------------------------------------------------------
# output postscript
#--------------------------------------------------------------------
    def output_postscript(self, filename):
        with open(filename, 'w') as f:
            bonds_output_postscript(self.bonds, f,
                    withforces=True, scaleforce=1.0/self.damp)

#--------------------------------------------------------------------
# run
#--------------------------------------------------------------------
    def run(self, Niter):
        bonds = self.bonds
        pairs = self.pairs
        groups = self.groups
        i0 = self.it

        # start iteration
        start = time.time()

        DAMP = self.damp
        filemode = 'w' if i0 == 0 else 'a'
        fdist = open('spadic10_revA.log', filemode)

        try:
         for i in range(i0, i0+Niter):
             # output various minimum distances TODO make it faster
             s = '%4i %6.1f %6.1f %6.1f %6.1f' % (i,
               min(abs(pair._dist_pboard()) for pair in pairs
                   if all(bond.group == 4 for bond in pair.bonds)),
               min(abs(pair._dist_pboard()) for pair in pairs
                   if all(bond.group == 6 for bond in pair.bonds)),
               min([abs(pair._dist_pboard_wire(0)[0])
                    for pair in pairs if abs(pair._dist_pboard_wire(0)[1]-0.5) < 0.5] +
                   [abs(pair._dist_pboard_wire(1)[0])
                    for pair in pairs if abs(pair._dist_pboard_wire(1)[1]-0.5) < 0.5]),
               min([abs(pair._dist_pchip_wire(0)[0])
                    for pair in pairs if abs(pair._dist_pchip_wire(0)[1]-0.5) < 0.5] +
                   [abs(pair._dist_pchip_wire(1)[0])
                    for pair in pairs if abs(pair._dist_pchip_wire(1)[1]-0.5) < 0.5])
               )
             print s; print >> fdist, s

             # add repulsive forces
             for pair in pairs:
                 pair.repulsion_pboard(damp=1.0*DAMP)
                 pair.repulsion_pboard_wire(damp=0.55*DAMP)
                 pair.repulsion_pchip_wire(damp=1.0*DAMP)

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

             # add random forces
             for bond in bonds:
                 ftot = sum(bond.forces)
                 frandmax = abs(ftot)*0.2
                 frandx = 2*frandmax*(random.random()-0.5)
                 frandy = 2*frandmax*(random.random()-0.5)
                 frand = Point2D(frandx, frandy)
                 bond.add_force(frand)

             # apply all forces
             for bond in bonds:
                 bond.apply_force()

        except KeyboardInterrupt:
            print 'Stopped by user.'
                
        finally:
         stop = time.time()
         T = stop-start
         print 'took %.2f seconds (%.1fms per iteration)' % (T, 1000*T/Niter)

         # save state
         self.bonds = bonds
         self.it = i+1

         # print various resulting minimum distances
         print min(abs(pair._dist_pboard()) for pair in pairs
                   if all(bond.group == 4 for bond in pair.bonds))
         print min(abs(pair._dist_pboard()) for pair in pairs
                   if all(bond.group == 6 for bond in pair.bonds))

         # print bond wire lengths for each group
         for group in groups:
             print list(sorted(int(round(bond.length))
                               for bond in bonds if bond.group == group))

         # write postscript and altium output files
         self.output_postscript('spadic10_revA.ps')

         with open('spadic10_revA.pas', 'w') as f:
             bonds_output_altium(bonds, 'SPADIC10_revA', f)


#====================================================================
# MAIN
#====================================================================
if __name__=='__main__':
    bg = BondGenerator('input_spadic10_revA.txt')
    bg.run(20)

