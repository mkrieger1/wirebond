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
    def __init__(self, bonds, min_dist_pboard=None):
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


#==============================================================================
# functions for creating bond pairs
#==============================================================================
def all_pairs(bonds):
    pairs = []
    for i in range(len(bonds)):
        for j in range(i+1, len(bonds)):
            pairs.append(BondPair(bonds[i], bonds[j]))
    return pairs


def neighbor_pairs(bonds, rings, center=None):
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
    # len(pairs) should now be len(bonds) * len(rings)
    return pairs


#==============================================================================
# functions for doing something with bond pairs
#==============================================================================
def mark_for_update(pairs):
    for pair in pairs:
        if any(b.forces for b in pair.bonds):
            pair._needs_update = True


#==============================================================================
# geometric calculations
#==============================================================================
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

