import math
from bondGeneratorLib import Point2D, Bond, BondPair

def valid_lines(filename):
    with open(filename) as f:
        for line_raw in f:
            line = line_raw.strip()
            if (not line.startswith('#') and len(line)):
                yield line

def read_bond_definition(filename):
    pos = Point2D(0, 0)
    incr = Point2D(0, 0)
    bonds = []
    padnumber = 0
    read_pads = False

    # process input file
    for line in valid_lines(filename):
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

    return (bonds, rings)

