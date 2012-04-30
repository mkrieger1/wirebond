# -*- coding: utf-8 -*-

#==============================================================================
# This is a clone/port of Peter Fischer's bond generator program
# from C++ to Python
#
# changelog:
#
# 2012/02/10 initial clone: same functionality without GUI,
#            aimed for compatibility with input files
#
# 2012/02/15 changed format of input file, can now handle arbitrary geometries
#            to do: update algorithm accordingly
#==============================================================================

import math
import numpy as np


#==============================================================================
# Point2D: 2-dimensional vector
#==============================================================================
class Point2D():
    def __init__(self, x=0, y=0):
        self.x = float(x)
        self.y = float(y)

    def __str__(self): # --> "str(p)" returns "p.__str__()"
        return "(%7.1f, %7.1f)" % (self.x, self.y)

    def __abs__(self):
        return math.sqrt(self.x**2 + self.y**2)

    def __add__(self, p): # overload "+" operator
        x = self.x + p.x
        y = self.y + p.y
        return Point2D(x, y)

    def __sub__(self, p): # overload "-" operator
        x = self.x - p.x
        y = self.y - p.y
        return Point2D(x, y)

    def __mul__(self, a): # overload "*" operator (Point2D * a)
        x = a*self.x
        y = a*self.y
        return Point2D(x, y)
    __rmul__ = __mul__ # reverse (a * Point2D)

    def __neg__(self): # "-Point2D"
        return -1*self

    def __div__(self, a): # overload "/" operator
        x = self.x/float(a)
        y = self.y/float(a)
        return Point2D(x, y)

    def dot(self, other):
        return self.x*other.x + self.y*other.y

    def rotate(self, phi):
        x = math.cos(phi)*self.x - math.sin(phi)*self.y
        y = math.sin(phi)*self.x + math.cos(phi)*self.y
        self.x = x
        self.y = y

    def rotated(self, phi):
        p = self
        p.rotate(phi)
        return Point2D(p.x, p.y)

    def polar_angle(self):
        return math.atan2(self.y, self.x) % (2*math.pi)

    def normalize(self):
        self.x /= abs(self)
        self.y /= abs(self)

    def normalized(self):
        return self/abs(self)


#==============================================================================
# Bond: wire between two endpoints (pchip --> pboard), forces move pboard
#==============================================================================
class Bond():
    def __init__(self, padnumber, pchip, length, angle, ring=0):
        self.padnumber = padnumber
        self.pchip = pchip
        self.length = float(length)
        self.angle = float(angle)
        self.calc_pboard()
        self.ring = int(ring)
        self.forces = []

    def __str__(self):
        return ("Bond #%i on ring %i, %s --> %s, "
                "angle = %5.1f°, length = %6.1f µm") % (self.padnumber,
                self.ring, str(self.pchip), str(self.pboard),
                self.angle*180/math.pi, self.length)

    def calc_pboard(self):
        wire = self.length*Point2D(1, 0).rotated(self.angle)
        self.pboard = self.pchip + wire

    def calc_length_angle(self):
        wire = self.pboard - self.pchip
        self.length = abs(wire)
        self.angle = wire.polar_angle()

    def set_pboard(self, pboard):
        self.pboard = pboard
        self.calc_length_angle()

    def set_length(self, length):
        self.length = length
        self.calc_pboard()

    def set_angle(self, angle):
        self.angle = angle
        self.calc_pboard()

    def stretch(self, dist):
        self.set_length(self.length + dist)

    def rotate(self, phi):
        self.set_angle((self.angle + phi) % (2*math.pi))

    def move(self, vector):
        self.set_pboard(self.pboard + vector)

    #def polar_offset(self):
    #    dphi = self.angle - self.pchip.polar_angle()
    #    return (dphi + math.pi) % (2*math.pi) - math.pi

    def add_force(self, force):
        self.forces.append(force)

    def apply_force(self):
        # start summation from Point2D(0, 0), not from int(0)
        force = sum(self.forces, Point2D(0,0))
        self.forces = []
        # add force to end point and rotate to that angle
        self.set_angle((self.pboard - self.pchip + force).polar_angle())


#==============================================================================
# BondPair: pair of two bonds
#==============================================================================
class BondPair():
    def __init__(self, bonds, min_dist_pboard=300):
        self.bonds = bonds
        self.min_dist_pboard = min_dist_pboard
        self._needs_update = True

    def dist_pboard(self):
        if self._needs_update:
            self._dist_pboard = abs(self.bonds[1].pboard -
                                    self.bonds[0].pboard)
            self._needs_update = False
        return self._dist_pboard

    def repulsion_pboard(self, damp=1.0):
        if self.dist_pboard() < self.min_dist_pboard:
            dir_ = (self.bonds[0].pboard - self.bonds[1].pboard).normalized()
            f = dir_ * (self.min_dist_pboard-self.dist_pboard()) * 0.5 * damp
            self.bonds[0].add_force(f)
            self.bonds[1].add_force(-f)

    #def update(self):
    #    self._dir_p1 = self.pair[1].pchip - self.pair[0].pchip
    #    self.dist_p1 = abs(self._dir_p1)
    #    self._dir_p2 = self.pair[1].pboard - self.pair[0].pboard
    #    self.dist_p2 = abs(self._dir_p2)
    #    self._needs_update = False

    #def repulsion_p1_l(self, damp=1.0):
    #    for p in [self.pair, list(reversed(self.pair))]:
    #        bond1 = p[0]
    #        bond2 = p[1]
    #        a = bond2.pchip-bond1.pchip
    #        b = bond1.pboard-bond1.pchip
    #        q = bond1.pchip + b.normed()*a.dot(b)/abs(b)
    #        print q
    #        F = q - bond2.pchip
    #        if abs(F) < 90:
    #            f = abs(q-bond1.pchip)/bond1.length*F.set_length(90-abs(F))
    #            bond1.add_force(damp*f)

        

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
def bonds_intersect(bond1, bond2):
    parallel = not ((bond2.phi - bond1.phi) % math.pi)
    if not parallel:
        l1   = bond1.l;    l2   = bond2.l
        phi1 = bond1.phi;  phi2 = bond2.phi
        x1   = bond1.p1.x; x2   = bond2.p1.x
        y1   = bond1.p1.y; y2   = bond2.p1.y

        # A.q = p
        A = np.array([[l1*math.cos(phi1), -l2*math.cos(phi2)],
                      [l1*math.sin(phi1), -l2*math.sin(phi2)]])
        p = np.array([[x2-x1],
                      [y2-y1]])

        # q = Ainv.p
        Ainv = np.linalg.inv(A)
        q = np.dot(Ainv, p)

        # ti: position of intersection on bond_i
        # ti = 0 --> p1 (pad on chip)
        # ti = 1 --> p2 (pad on PCB)
        (t1, t2) = (q[0][0], q[1][0])

    else: # parallel
        (t1, t2) = (None, None)

    return (t1, t2)



#------------------------------------------------------------------------------
# Step 1 - Read chip pad geometry from file and create a "Bond" for each pad
#------------------------------------------------------------------------------
def get_valid_lines(filename):
    with open(filename) as f:
        for line_raw in f:
            line = line_raw.strip()
            if (not line.startswith('#') and len(line)):
                yield line

def read_chip_pad_definitions(filename):
    pos = Point2D(0, 0)
    incr = Point2D(0, 0)
    bonds = []
    padnumber = 0
    read_pads = False

    # process input file
    for line in get_valid_lines(filename):
        # get rings
        if line.startswith('RINGS'):
            rings = map(int, line.split()[1:])

        # get min./max. angle
        elif line.startswith('ANGLES'):
            angles_int = map(int, line.split()[1:])
            angles_fixed = not(abs(angles_int[1]-angles_int[0]) == 360)
            angles = [x*math.pi/180 for x in angles_int]

        # get minimum pad distance
        elif line.startswith('D_MIN'):
            d_min = int(line.split()[1])

        # move position
        elif line.startswith('MOVE'):
            [dx, dy] = map(int, line.split()[1:])
            pos.x += dx
            pos.y += dy

        # change position auto-increment
        elif line.startswith('INCR'):
            [dx, dy] = map(int, line.split()[1:])
            incr.x = dx
            incr.y = dy

        # go into read_pad mode
        elif line.startswith('BEGIN PADS'):
            read_pads = True

        # leave read_pad mode
        elif line.startswith('END PADS'):
            read_pads = False
            pos -= incr # undo position change after last pad

        # read pads if in read_pad mode and
        # create a bond for each pad with ring > -1
        elif read_pads:
            pads = map(int, line.split())
            for pad in pads:
                ring = pad-1
                if ring > -1: # ring == -1 is an empty pad
                    pchip = Point2D(pos.x, pos.y)
                    length = rings[ring]
                    bonds.append(Bond(padnumber, pchip, length, 0, ring))
                pos += incr
                padnumber += 1

    # distribute bonds evenly between min. and max. angle
    step_phi = (angles[1]-angles[0]) / (len(bonds)-1)
    phi = angles[0]
    for b in bonds:
        b.set_angle(phi)
        phi += step_phi
    #if angles_fixed:
    #    bonds[ 0].angle_fixed = True
    #    bonds[-1].angle_fixed = True
    #    print "fixed bond %s to %.1f degrees" % (bonds[0].name, bonds[0].phi*180/math.pi)
    #    print "fixed bond %s to %.1f degrees" % (bonds[-1].name, bonds[-1].phi*180/math.pi)
    #else:
    #    print "angles not fixed"

    return (rings, bonds)


#------------------------------------------------------------------------------
# Step 2 - Iteratively shift bonds around
#------------------------------------------------------------------------------
def create_all_bondpairs(bonds):
    pairs = []
    for i in range(len(bonds)):
        for j in range(i+1, len(bonds)):
            min_dist_pboard = 90 if (bonds[i].ring == bonds[j].ring and
                                     bonds[i].ring in [0, 1, 2]) else 300
            pairs.append(BondPair(bonds[i], bonds[j], min_dist_pboard))
    return pairs

def neighbor_pairs(rings, bonds, center=None):
    pairs = []
    # sort all bonds clockwise with respect to their chip pad position use
    # center of gravity or custom center point (custom is useful if the bond
    # pads are not distributed evenly around the whole chip, e.g. only one
    # edge)
    if center is None:
        center = sum((b.pchip for b in bonds), Point2D(0, 0))/len(bonds)
    bonds_cw = list(sorted(bonds,
                    key=lambda b: (b.pchip-center).polar_angle()))
    # for each bond, collect the next other bond on each ring (do this only
    # in clockwise direction to avoid duplicate pairs)
    for (i, bond) in enumerate(bonds_cw):
        neighbors = [None for ring in rings]
        otherbonds = iter(bonds_cw[i+1:] + bonds_cw[:i])
                        # i+1, i+2, ..., 0, 1, ..., i-1
        while not all(neighbors): # bool(None) == False, bool(any bond) == True
            otherbond = otherbonds.next()
            ring = otherbond.ring
            if neighbors[ring] is None:
                neighbors[ring] = otherbond
        for neighbor in neighbors:
            pairs.append(BondPair([bond, neighbor]))
    # len(pairs) == len(rings) * len(bonds)
    return pairs



def mark_processed_pairs(pairs):
    for p in pairs:
        if any(len(b.forces) for b in p.pair):
            p._needs_update = True

def process_pairs_p2(pairs_p2):
    # sort pairs_in_range by their p2 distance
    pairs_p2.sort() # only needed for displaying min. distance
    p = pairs_p2[0]
    print "min. p2 distance:", str(p), p.dist_p2,

    # add repulsive forces if necessary
    for p in pairs_p2:
        if p.dist_p2 < p.min_dist_pboard:
            p.repulsion_p2()

def process_all_bonds(bonds, pairs, pairs_p2):
    process_pairs_p2(pairs_p2)

    # buggy
    #for p in pairs:
    #    p.repulsion_p1_l()

    mark_processed_pairs(pairs)

    for b in bonds: b.apply_force()

    updated = 0
    for p in pairs_p2:
        if p._needs_update:
            p.update()
            updated += 1
    print updated, "pairs updated"

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

