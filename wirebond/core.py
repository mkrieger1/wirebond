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
        if p == 0: # Point2D + 0
            return self
        else:
            x = self.x + p.x
            y = self.y + p.y
            return Point2D(x, y)
    __radd__ = __add__ # 0 + Point2D()

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
    def __init__(self, padnumber, net, pchip, length, angle, ring,
                       rectangular, chip=None):
        self.padnumber = padnumber
        self.net = net
        self.pchip = pchip
        self.length = float(length)
        self.angle = float(angle)
        self.ring = ring
        self.rectangular = rectangular

        if self.rectangular:
            self.ring_radius = self.length
            self.chip = chip

        self.pad_rotation = self.angle
        self.calc_pboard()
        self.forces = []

    def __str__(self):
        return ("Bond #%i on ring %i, %s --> %s, "
                "angle = %5.1f°, length = %6.1f µm, net = %s") % (self.padnumber,
                self.ring, str(self.pchip), str(self.pboard),
                self.angle*180/math.pi, self.length, self.net)

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
        force = sum(self.forces) # works because 0 + Point2D is implemented
        self.forces = []
        # add force to end point and rotate to that angle
        self.set_angle((self.pboard - self.pchip + force).polar_angle())
        self.pad_rotation = self.angle
        # if board-side pad sits on a rectangle around the chip instead of a
        # circle around the chip-side pad, adjust the length of the bond
        if self.rectangular:
            self.set_to_rectangle(rounded_corners=True)

    def set_to_rectangle(self, rounded_corners):
        R = 0 # experimental
        [x0, y0, x1, y1] = self.chip
        radius = self.ring_radius + R
        wire = [self.pchip, self.pboard]
        quadrant = int(self.angle*2/math.pi)

        right = quadrant in [0, 3]
        top   = quadrant in [0, 1]
        (x, sgnx) = (x1, 1) if right else (x0, -1)
        (y, sgny) = (y1, 1) if top   else (y0, -1)
        x -= sgnx * R
        y -= sgny * R
        corner = Point2D(x, y)

        t = min(t for t in [intersect_line_x(wire, x + sgnx*radius)[1],
                            intersect_line_y(wire, y + sgny*radius)[1]]
                if t > 0)

        self.set_length(t*self.length)
        self.pad_rotation = int((self.angle+math.pi/4)*2/math.pi)*math.pi/2

        if rounded_corners:
            if (sgnx * (self.pboard.x - corner.x) > 0 and
                sgny * (self.pboard.y - corner.y) > 0):
             wire = [self.pchip-corner, self.pboard-corner]
             (q, t) = intersect_line_circle(wire, radius, quadrant)
             self.set_length(t*self.length)
             self.pad_rotation = q.polar_angle()


#==============================================================================
# BondPair: pair of two bonds
#==============================================================================
class BondPair():
    def __init__(self, bonds, min_dist_pboard=None,
                              min_dist_pboard_wire=None,
                              min_dist_pchip_wire=None):
        self.bonds = bonds
        self.set_min_dist_pboard(min_dist_pboard)
        self.set_min_dist_pboard_wire(min_dist_pboard_wire)
        self.set_min_dist_pchip_wire(min_dist_pchip_wire)

    def _bonds_perm(self, index):
        return [self.bonds[index], self.bonds[1-index]]

    def _dist_pboard(self):
        [bond, otherbond] = self._bonds_perm(0)
        return bond.pboard-otherbond.pboard

    def _dist_pboard_wire(self, index):
        [bond, otherbond] = self._bonds_perm(index)
        return dist_point_line(otherbond.pboard, [bond.pchip, bond.pboard])

    def _dist_pchip_wire(self, index):
        [bond, otherbond] = self._bonds_perm(index)
        return dist_point_line(otherbond.pchip, [bond.pchip, bond.pboard])

    def set_min_dist_pboard(self, value):
        self._min_dist_pboard = value

    def set_min_dist_pboard_wire(self, value):
        self._min_dist_pboard_wire = value

    def set_min_dist_pchip_wire(self, value):
        self._min_dist_pchip_wire = value

    def repulsion_pboard(self, damp=1.0):
        [bond, otherbond] = self._bonds_perm(0)
        dist = self._dist_pboard()
        dist_violation = self._min_dist_pboard - abs(dist)
        if dist_violation > 0:
            f = damp * dist_violation * dist.normalized()
            bond.add_force(0.5*f)
            otherbond.add_force(-0.5*f)

    def repulsion_pboard_wire(self, damp=1.0):
        for i in range(2):
            [bond, otherbond] = self._bonds_perm(i)
            (dist, t) = self._dist_pboard_wire(i)
            dist_violation = self._min_dist_pboard_wire - abs(dist)
            if t > 0 and t < 1 and dist_violation > 0:
                f = damp * dist_violation * dist.normalized()
                bond.add_force(0.5*f)
                otherbond.add_force(-0.5*f)

    def repulsion_pchip_wire(self, damp=1.0):
        for i in range(2):
            [bond, otherbond] = self._bonds_perm(i)
            (dist, t) = self._dist_pchip_wire(i)
            dist_violation = self._min_dist_pchip_wire - abs(dist)
            if t > 0 and t < 1 and dist_violation > 0:
                f = damp * dist_violation * dist.normalized()
                bond.add_force(f)


#==============================================================================
# functions for creating bond pairs
#==============================================================================
def all_pairs(bonds):
    pairs = []
    for i in range(len(bonds)):
        for j in range(i+1, len(bonds)):
            pairs.append(BondPair([bonds[i], bonds[j]]))
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
# geometric calculations
#==============================================================================

#------------------------------------------------------------------------------
# distance point <--> line
#------------------------------------------------------------------------------
def dist_point_line(point, line):
    '''Calculate the distance between a point and a line.

    The line is defined by two Point2D instances given in a list [p1, p2],
    the point is given by a Point2D instance p.
    
    Returns a tuple (r, t) where b = p1 + t * (p2-p1) is the point on the
    line with the shortest distance to p, and r is the vector from p to b.
    The distance between the point and the line is then abs(r).'''
    a = point-line[0]
    b = (line[1]-line[0]).normalized()
    t = a.dot(b)
    return (line[0] + b*t - point, t/abs(line[1]-line[0]))

#------------------------------------------------------------------------------
# intersection line <--> line
#------------------------------------------------------------------------------
def _intersect_line_xy(line, xy, mode):
    p = line[0]
    a = line[1]-line[0]
    [P, A] = {
      'x': [p.x, a.x],
      'y': [p.y, a.y]
      }[mode]
    if A == 0: # line parallel to x or y
        if P == xy: # number of solutions: infinite -> choose one
            t = 1.0
        else:       # number of solutions: 0
            t = None
    else:
        t = (xy-P)/A
    q = p + t*a if t is not None else None
    return (q, t)

def intersect_line_x(line, x):
    '''Find the point on a line with a given x coordinate.
    
    The line is defined by two Point2D instances given in a list [p1, p2].
    
    Returns a tuple (q, t) where q = p1 + t*(p2-p1) is the intersection
    point. If the line itself is the vertical line defined by x, from the
    infinite number of solutions the one where t=1 is chosen. If the line is
    a vertical line with a different x coordinate, q and t are None.'''
    return _intersect_line_xy(line, x, 'x')

def intersect_line_y(line, y):
    '''Find the point on a line with a given y coordinate.

    The line is defined by two Point2D instances given in a list [p1, p2].
    
    Returns a tuple (q, t) where q = p1 + t*(p2-p1) is the intersection
    point. If the line itself is the horizontal line defined by y, from the
    infinite number of solutions the one where t=1 is chosen. If the line is
    a horizontal line with a different y coordinate, q and t are None.'''
    return _intersect_line_xy(line, y, 'y')

#------------------------------------------------------------------------------
# intersection line <--> circle
#------------------------------------------------------------------------------
def intersect_line_circle(line, radius, quadrant):
    '''Calculate the intersection point of a line and a quarter circle
    centered at the origin.
    
    source: http://mathworld.wolfram.com/Circle-LineIntersection.html
    
    The line is defined by two Point2D instances given in a list [p1, p2].
    The quadrant is given as the number 0, 1, 2, or 3 where 0 is top right,
    1 is top left, 2 is bottom left and 3 is bottom right,
    or quadrant = int(phi*2/pi).
    
    Returns a tuple (q, t) where q = p1 + t*(p2-p1) is the intersection
    point. If there is no intersection point, q and t are None.'''

    [p1, p2] = line
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    sgn = -1 if dy < 0 else 1
    l = abs(p2-p1)
    if l == 0:
        raise ValueError('line is undefined, need two different points')
    D = p1.x*p2.y - p1.y*p2.x
    DD = (radius*l)**2 - D**2

    if DD < 0:
        result = []
    else:
        result = [q
          for q in [Point2D(( D*dy + sgn*dx  * math.sqrt(DD)) / l**2,
                            (-D*dx + abs(dy) * math.sqrt(DD)) / l**2),
                    Point2D(( D*dy - sgn*dx  * math.sqrt(DD)) / l**2,
                            (-D*dx - abs(dy) * math.sqrt(DD)) / l**2)]
          if int(q.polar_angle()*2/math.pi) == quadrant]

    if result:
        q = result[0]
        t = abs(q-p1)/l
    else:
        q = None
        t = None
    return (q, t)
        
        


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

