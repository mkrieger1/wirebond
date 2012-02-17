from bondGeneratorLib import *

# read input file
(rings, angles, d_min, bonds) = \
  read_chip_pad_definitions('input/spadic10_pads.txt')

# testing bond manipulating algorithms
pitch = 95
(bonds, d_av, d_min) = \
  iterate_bonds(rings, angles, pitch, bonds, 100, shiftInterposer=False)
#bonds = clean_up(bonds, rings, d_av)


from bondOutputPostscript import *
with open('bondspstest.ps', 'w') as f:
    bonds_output_postscript(bonds, f)

from bondOutputAltium import *
with open('spadic10.pas', 'w') as f:
    bonds_output_altium(bonds, 'spadic10', f)

