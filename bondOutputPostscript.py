#------------------------------------------------------------------------------
# Beginning of Postscript file: header & subroutine definitions
#------------------------------------------------------------------------------
def postscript_header():
    return """%!PS-Adobe-3.0 EPSF-3.0
%%BoundingBox: 0 0 595 842

/COL_PAD    {0.5 0.8   0 setrgbcolor} def
/COL_BOND   {0   dup dup setrgbcolor} def

/PADSml_W 50 def
/PADSml_L 90 def
/PADLrg_W 150 def
/PADLrg_L 300 def

/Pad_Small {      %% on stack: rot x y
  gsave translate rotate
  PADSml_L -2 div PADSml_W -2 div PADSml_L PADSml_W rectfill
  grestore
} def

/Pad_Large {      %% on stack: rot x y
  gsave translate rotate
  PADLrg_L -2 div PADLrg_W -2 div PADLrg_L PADLrg_W rectfill
  grestore
} def

/Bond {           %% on stack: x1 y1 x2 y2 phi
  5 dict begin
  /phi exch def
  /y2  exch def
  /x2  exch def
  /y1  exch def
  /x1  exch def
  gsave

  COL_PAD
  phi x1 y1 Pad_Small
  phi x2 y2 Pad_Large

  COL_BOND 10 setlinewidth
  x1 y1 moveto x2 y2 lineto stroke

  grestore
  end
} def

/DrawBonds {"""


#------------------------------------------------------------------------------
# Draw a bond
#------------------------------------------------------------------------------
def postscript_bond(bond):
    return "%.1f %.1f %.1f %.1f %.1f Bond" % (bond.phi, bond.p1.x, bond.p1.y,
                                                        bond.p2.x, bond.p2.y)


#------------------------------------------------------------------------------
# End of Postscript file
#------------------------------------------------------------------------------
def mm2pt(x):
    return x/25.4*72

def postscript_footer(scale=10):
    return """} def

%.6f dup scale
200 300 translate
DrawBonds

showpage 
""" % mm2pt(1.0*scale) # e.g. scale=10 --> 100 um real = 1 mm on paper


#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
def bonds_output_postscript(bonds):
    print postscript_header()
    for bond in bonds:
        print postscript_bond(bond)
    print postscript_footer(scale=10)

