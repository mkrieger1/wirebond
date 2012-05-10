import math

#------------------------------------------------------------------------------
# Beginning of Postscript file: header & subroutine definitions
#------------------------------------------------------------------------------
def postscript_header():
    return """%!PS-Adobe-3.0 EPSF-3.0
%%BoundingBox: 0 0 595 842

%------------------------------------------------------------------------------
% customize here
%------------------------------------------------------------------------------
/COL_PAD    {0.9 0.8 0.5 setrgbcolor} def
/COL_BOND   {0   0   0   setrgbcolor} def
/COL_FORCE  {0.8 0.1 0.1 setrgbcolor} def

/PADSml_W  50 def
/PADSml_L  90 def
/PADLrg_W 150 def
/PADLrg_L 150 def

/Pad_Small { % on stack: rot x y
  gsave translate rotate
  PADSml_L -2 div PADSml_W -2 div PADSml_L PADSml_W rectfill
  grestore
} def

/Pad_Large { % on stack: rot x y
  gsave translate rotate
  PADLrg_L -2 div PADLrg_W -2 div PADLrg_L PADLrg_W rectfill
  grestore
} def

/Bond { % on stack: x1 y1 x2 y2 phi
  gsave 5 dict begin
  /phi exch def
  /y2  exch def
  /x2  exch def
  /y1  exch def
  /x1  exch def

  COL_PAD
  phi x1 y1 Pad_Small
  phi x2 y2 Pad_Large

  COL_BOND 10 setlinewidth
  x1 y1 moveto x2 y2 lineto stroke

  end grestore
} def

/Force { % on stack: x1 y1 x2 y2
  gsave 4 dict begin
  /y2 exch def
  /x2  exch def
  /y1  exch def
  /x1  exch def

  COL_FORCE 10 setlinewidth
  x1 y1 moveto x2 y2 lineto stroke

  end grestore
} def

%------------------------------------------------------------------------------
% scaling helper functions
%------------------------------------------------------------------------------
/per_mm {
  % x per_mm --> x units = 1 mm on paper
  72 25.4 div exch div dup scale
} def

/mm {
  % x mm --> 72/25.4*x pt
  72 25.4 div mul
} def


%------------------------------------------------------------------------------
% list of all bonds
%------------------------------------------------------------------------------
/DrawBonds { gsave"""



#------------------------------------------------------------------------------
# Draw a bond
#------------------------------------------------------------------------------
def postscript_bond(bond):
    return "%7.1f %7.1f %7.1f %7.1f %6.1f Bond" % (bond.pchip.x, bond.pchip.y,
                                                   bond.pboard.x, bond.pboard.y,
                                                   bond.angle/math.pi*180)


#------------------------------------------------------------------------------
# Draw a force line
#------------------------------------------------------------------------------
def postscript_force(bond, force, scaleforce=1.0):
    return "%7.1f %7.1f %7.1f %7.1f Force" % (bond.pboard.x, bond.pboard.y,
                                              bond.pboard.x + scaleforce*force.x,
                                              bond.pboard.y + scaleforce*force.y)


#------------------------------------------------------------------------------
# End of Postscript file
#------------------------------------------------------------------------------
def postscript_footer(per_mm=100):
    return """grestore } def

105 mm 120 mm translate

gsave
%.1f per_mm DrawBonds
grestore

gsave
1 setlinewidth
-5 mm 0 mm moveto 5 mm 0 mm lineto stroke
0 mm -5 mm moveto 0 mm 5 mm lineto stroke
grestore

showpage 
""" % per_mm


#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
def bonds_output_postscript(bonds, f=None, withforces=False, scaleforce=1.0):
    print >> f, postscript_header()
    for bond in bonds:
        print >> f, postscript_bond(bond)
    if withforces:
        for bond in bonds:
            for force in bond.forces:
                print >> f, postscript_force(bond, force, scaleforce)
    print >> f, postscript_footer(per_mm=100)

