#!/usr/bin/python
# -*- coding: utf-8 -*-

#------------------------------------------------------------------------------
# This is a clone/port of Peter Fischer's bond generator program
# from C++ to Python
#
# changelog:
#
# 2012/02/10 initial clone: same functionality without GUI,
#            aimed for compatibility with input files
#------------------------------------------------------------------------------

import math

#------------------------------------------------------------------------------
# Classes
#------------------------------------------------------------------------------

class Point2D():
    def __init__(self, x=0, y=0):
        self.x = float(x)
        self.y = float(y)

    def __str__(self): # --> "str(p)" returns "p.__str__()"
        return "(%7.1f, %7.1f)" % (self.x, self.y)

class Bond():
    def __init__(self, name, p1, p2, phi=0, l=0, row=0):
        self.name = name
        self.p1 = p1
        self.p2 = p2
        self.phi = float(phi)
        self.l = float(l)
        self.row = int(row)

    def calc_p2_from_lphi(self):
        self.p2.x = self.p1.x + self.l * math.cos(self.phi)
        self.p2.y = self.p1.y + self.l * math.sin(self.phi)

    def calc_lphi_from_p2(self):
        dx = self.p2.x-self.p1.x
        dy = self.p2.y-self.p1.y
        self.phi = math.atan2(dy, dx)
        self.l = math.sqrt(dx**2 + dy**2)

    def rotate_angle(self, dphi):
        self.phi += dphi
        self.calc_p2_from_lphi()

    def rotate_dist(self, d):
        self.rotate_angle(d/self.l)

    def __str__(self):
        return ("Bond \"%s\" on row %i, %s --> %s, "
                "phi = %5.1f°, l = %6.1f µm") % (self.name, self.row,
                str(self.p1), str(self.p2), self.phi*180/math.pi, self.l)

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
def distance(bond, point2d):
    b = bond
    a = point2d
    dirx = math.cos(b.phi)
    diry = math.sin(b.phi)
    lambd = (b.p2.x - a.x) * dirx + (b.p2.y - a.y) * diry
    dx = a.x - b.p2.x + lambd * dirx
    dy = a.y - b.p2.y + lambd * diry
    return math.sqrt(dx**2+dy**2) if (lambd > 0) else 0

def distance2(bond, point2d):
    b = bond
    a = point2d
    dirx = math.cos(b.phi)
    diry = math.sin(b.phi)
    dx = (b.p2.x - a.x)
    dy = (b.p2.y - a.y)
    return abs(dx*diry - dy*dirx)

def pointDistance(bond1, bond2):
    return max(distance(bond1, bond2.p2), distance(bond2, bond1.p2))


#------------------------------------------------------------------------------
# Step 1 - Read chip pad geometry from file and create a "Bond" for each pad
#------------------------------------------------------------------------------
def read_chip_pad_definitions(filename):
    f = open(filename)
    
    # read number of rows and radii
    ringline = map(int, f.readline().split())
    Nrings = ringline[0] # TO DO: change in input file format
                         # do not specify number of rings, only their radii
    rings = ringline[1:] # --> Nrings = len(rings)

    # read angles
    angles = map(int, f.readline().split())
    PHI_L = (angles[0]+90) * math.pi/180
    PHI_R = (90-angles[1]) * math.pi/180

    # read pitch
    pitch = int(f.readline()) # um

    # read chip pads
    pads = ' '.join(f.readlines()[:-1]).split() # last line contains '-1' as EOF marker
    bonds = []
    xpos = 0

    # create bonds
    for (padnumber, pad) in enumerate(pads):
        row = int(pad)-1
        if row >= 0: # row == -1 is an empty pad
            p1 = Point2D(xpos, 0)
            p2 = Point2D(0, 0)
            l = rings[row]
            bonds.append(Bond(str(padnumber), p1, p2, 0, l, row))
        xpos += pitch

    # distribute evenly between min. and max. angle
    step_phi = (PHI_R - PHI_L) / (len(bonds)-1)
    phi = PHI_L
    for b in bonds:
        b.phi = phi
        b.calc_p2_from_lphi()
        phi += step_phi

    f.close()

    return (rings, angles, pitch, bonds)


#------------------------------------------------------------------------------
# Step 2 - Iteratively shift bonds around
#------------------------------------------------------------------------------
def iterate_bonds(rings, angles, pitch, bonds, NITER, shiftInterposer=False):
    DMIN = pitch
    N = len(bonds)

    for it in range(NITER):
        if shiftInterposer:
            for i in range(1, N):
                pass
#               float d12 = PointDistance(BB[i-1], BB[i  ]);  // shift on interposer
#               float d23 = PointDistance(BB[i  ], BB[i+1]);
#               float dav = (d12 + d23) / 2.0;
#               float dx = PI[i].x - PI[i-1].x;
#               PI[i].x = PI[i-1].x + dx * dav / d12;
        
        d_av = 0
        d_min = -1
        npair = 0

        for i in range(1, N):
            iprev = i-1
            if (bonds[iprev].row == bonds[i].row):
                dist = pointDistance(bonds[i-1], bonds[i])
                d_av += dist
                d_min = dist if (d_min == -1) else min(d_min, dist)
                npair += 1
        if npair > 0:
            d_av /= npair

        for i in range(1, N-1):
            inext = i+1
            iprev = i-1
            prev_is_same_row = (bonds[iprev].row == bonds[i].row)
            next_is_same_row = (bonds[inext].row == bonds[i].row)
            d12 = pointDistance(bonds[iprev], bonds[i])
            d23 = pointDistance(bonds[i], bonds[inext])
            d13 = d12 + d23

            if (prev_is_same_row and next_is_same_row):
                dav = d13 / 2.0
                shift = d12 - dav
                bonds[i].rotate_dist(shift)

            else:
                d12_target = d_av if prev_is_same_row else DMIN
                d23_target = d_av if next_is_same_row else DMIN
                d13_target = d12_target + d23_target

                shift1 = 0 if (iprev == 0) else (d13 - d13_target) / 2.0
                shift3 = -(d13 - d13_target - shift1)
                shift2 = shift1 + d12_target - d12

                if (inext == N-1):
                    shift1 -= shift3
                    shift2 -= shift3
                    shift3 -= shift3 # 0

                bonds[iprev].rotate_dist(-shift1)
                bonds[i    ].rotate_dist(-shift2)
                bonds[inext].rotate_dist(-shift3)

    print "PCB Bond average: ", d_av
    print "PCB Bond min    : ", d_min

    return (bonds, d_av, d_min)


#------------------------------------------------------------------------------
# Step 3 - Clean Up (better do it right in the first place)
#------------------------------------------------------------------------------
def clean_up(bonds, rings, d_av):
    N = len(bonds)
    OUTERROW = len(rings)-1
    i_last = -1
    iproblem = -1

    for i in range(N-1):
        if (i_last == -1):
            if (bonds[i].row == OUTERROW):
                i_last = i
        else:
            if (bonds[i].row == OUTERROW):
                dx = bonds[i].p2.x - bonds[i_last].p2.x
                dy = bonds[i].p2.y - bonds[i_last].p2.y
                d_bad = math.sqrt(dx**2 + dy**2)
                if (dx < 0):
                    iproblem = i
                    break
                i_last = i

    if (iproblem == -1):
        print "Nothing to clean up"
    else:
        print "Cleaning up index ", iproblem
        print "d_av = ", d_av
        d_bad *= -1
        print "Bad distance = ", d_bad
        shift = (d_av - d_bad) / 2.0
        angle = shift / rings[OUTERROW]

        for i in range(iproblem):
            bonds[i].rotate_angle(angle * i / (iproblem-1))
        for i in range(iproblem, N):
            bonds[i].rotate_angle(-angle * (N-1-i) / (N-1-iproblem))

    return bonds


#------------------------------------------------------------------------------
# Step 4 - Create XY File
#------------------------------------------------------------------------------

# do this in separate module



#------------------------------------------------------------------------------
# Step 5 - Create Postscript File
#------------------------------------------------------------------------------

# do this in separate module

        

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
if __name__=='__main__':
    (rings, angles, pitch, bonds) = \
      read_chip_pad_definitions('input/Pattern_Top.txt')
    (bonds, d_av, d_min) = \
      iterate_bonds(rings, angles, pitch, bonds, 100, shiftInterposer=False)
    bonds = clean_up(bonds, rings, d_av)

    from bondOutputPostscript import *
    with open('bondspstest.ps', 'w') as f:
        bonds_output_postscript(bonds, f)

