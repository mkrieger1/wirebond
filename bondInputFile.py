import math
from bondGeneratorLib import Point2D, Bond, BondPair

def valid_lines(filename):
    with open(filename) as f:
        for line_raw in f:
            line = line_raw.strip()
            if (not line.startswith('#') and len(line)):
                yield line

def read_bond_definition(filename):
    rings = {}
    pos = Point2D(0, 0)
    incr = Point2D(0, 0)
    bonds = []
    padnumber = 0
    read_pads = False

    # process input file
    for line in valid_lines(filename):
        # get chip bounding box
        if line.startswith('CHIP'):
            chip = map(int, line.split()[1:])

        # get rings
        if line.startswith('RING'):
            item = line.split()
            ring = int(item[1])
            radius = int(item[2])
            rectangular = 'RECT' in item
            rings[ring] = (radius, rectangular)

        # get 'do not bond' nets
        if line.startswith('NO_BOND'):
            no_bond = line.split()[1:]

        # get min./max. angle
        #elif line.startswith('ANGLES'):
        #    angles_int = map(int, line.split()[1:])
        #    angles_fixed = not(abs(angles_int[1]-angles_int[0]) == 360)
        #    angles = [x*math.pi/180 for x in angles_int]

        # get minimum pad distance
        #elif line.startswith('D_MIN'):
        #    d_min = int(line.split()[1])

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
            item = line.split()
            padnumber = int(item[0])
            net = item[1]
            if net not in no_bond:
                ring = int(item[2])
                (length, rectangular) = rings[ring]
                pchip = Point2D(pos.x, pos.y)
                bond = Bond(padnumber, net, pchip, length, 0, ring,
                            rectangular, chip)
                bonds.append(bond)
            pos += incr

            #pads = map(int, line.split())
            #for pad in pads:
            #    ring = pad-1
            #    if ring > -1: # ring == -1 is an empty pad
            #        pchip = Point2D(pos.x, pos.y)
            #        length = rings[ring]
            #        bonds.append(Bond(padnumber, pchip, length, 0, ring))


    # distribute bonds evenly between min. and max. angle
    step_phi = 2*math.pi / (len(bonds)-1)
    phi = 3*math.pi/4
    #step_phi = (angles[1]-angles[0]) / (len(bonds)-1)
    #phi = angles[0]
    for b in bonds:
        b.set_angle(phi)
        b.apply_force() # in case rectangle is set, this sets the correct length
        phi += step_phi
    #if angles_fixed:
    #    bonds[ 0].angle_fixed = True
    #    bonds[-1].angle_fixed = True
    #    print "fixed bond %s to %.1f degrees" % (bonds[0].name, bonds[0].phi*180/math.pi)
    #    print "fixed bond %s to %.1f degrees" % (bonds[-1].name, bonds[-1].phi*180/math.pi)
    #else:
    #    print "angles not fixed"

    return (bonds, rings)

