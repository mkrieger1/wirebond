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
    rings = map(int, f.readline().split())
    Nrow = rings[0] # TO DO: change in input file format
                    # do not specify number of rings, only their radii
    Ring = rings[1:] # --> Nrow = len(Ring)

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
    phi = PHI_L
    step_phi = (PHI_R - PHI_L) / (len(pads)-1)

    # create bonds
    for (padnumber, pad) in enumerate(pads):
        p1 = Point2D(xpos, 0)
        p2 = Point2D(0, 0)
        row = int(pad)-1
        l = Ring[row]
        b = Bond(str(padnumber), p1, p2, phi, l, row)
        b.calc_p2_from_lphi()
        bonds.append(b)

        phi += step_phi
        xpos += pitch

    f.close()

    return (bonds, Ring)


#------------------------------------------------------------------------------
# Step 2 - ?
#------------------------------------------------------------------------------
def iterate_bonds(bonds, NITER, shiftInterposer=False):
    for it in range(NITER):
        if shiftInterposer:
            for i in range(1, len(bonds)):
                pass
#               float d12 = PointDistance(BB[i-1], BB[i  ]);  // shift on interposer
#               float d23 = PointDistance(BB[i  ], BB[i+1]);
#               float dav = (d12 + d23) / 2.0;
#               float dx = PI[i].x - PI[i-1].x;
#               PI[i].x = PI[i-1].x + dx * dav / d12;
        
        d_av = 0
        d_min = -1
        npair = 0

        for i in range(1, len(bonds)):
            if (bonds[i-1].row == bonds[i].row):
                dist = pointDistance(bonds[i-1], bonds[i])
                d_av += dist
                d_min = dist if (d_min == -1) else min(d_min, dist)
                npair += 1
        d_av /= npair
        

#//benoetigt! (Schritt 2)
#void __fastcall TForm1::ButtonIterateClick(TObject *Sender)
#{
#  Memo1->Clear();
#  const bool ShiftInterposer = CB_ShiftInterposer->Checked;
#  float d_min;
#
#  for (int it = 0; it < NITER->Text.ToInt(); it++) {
#
#    // Use pairs of distances and shift the middle point
#
#    if (ShiftInterposer) {
#      for (unsigned int i = 1; i < N-1; i++) {        // shift point i
#/*
#        float d12 = PointDistance(BB[i-1], BB[i  ]);  // shift on interposer
#        float d23 = PointDistance(BB[i  ], BB[i+1]);
#        float dav = (d12 + d23) / 2.0;
#        float dx = PI[i].x - PI[i-1].x;
#        PI[i].x = PI[i-1].x + dx * dav / d12;
#*/
#      }
#    }
#
#    // Calculate the average distance of adjacent bonds on board
#    d_av      =  0;
#    d_min     = -1;
#    int npair = 0;
#    for (unsigned int i = 1; i < N; i++) {
#      unsigned int iprev = i-1;
#      if (BB[iprev].row == BB[i].row) {
#        float dist = PointDistance(BB[iprev], BB[i    ]);
#        d_av += dist;
#        d_min = (d_min == -1) ? dist : std::min(d_min, dist);
#        npair++;
#      }
#    }
#    d_av /= npair;
#//    Memo1->Lines->Add("d_av: " + (AnsiString)d_av);
#
#    const float DMIN = pitch;
#
#    for (unsigned int i = 1; i < N-1; i++) {         // shift point i on board
#      unsigned int inext = i+1;
#      unsigned int iprev = i-1;
#      bool  prev_is_same_row = BB[iprev].row == BB[i].row;
#      bool  next_is_same_row = BB[i].row == BB[inext].row;
#      float d12 = PointDistance(BB[iprev], BB[i    ]);
#      float d23 = PointDistance(BB[i    ], BB[inext]);
#      float d13 = d12 + d23;
#
#      if (prev_is_same_row && next_is_same_row) {      // L+R same row
#        float dav = d13 / 2.0;                         // adjust middle point
#        float shift = d12 - dav;                       // linear shift
#        BB[i].Rotate_dist(shift);
#      } else {
#//        1        2                  3
#//        X        X                  X
#//               X   X        X              <-> d13_target
#//        |shift1|
#//        |           d13             |
#//               | d13_target |
#//                            |shift3 |
#
#        float d12_target = prev_is_same_row ? d_av : DMIN;
#        float d23_target = next_is_same_row ? d_av : DMIN;
#
#        float d13_target = d12_target + d23_target;
#        float shift1 = (iprev == 0) ? 0 : (d13 - d13_target) / 2.0;
#        float shift3 = - (d13 - d13_target - shift1);
#        float shift2 = shift1 + d12_target - d12;
#
#        if (inext == N-1) {
#          shift1 -= shift3;
#          shift2 -= shift3;
#          shift3 -= shift3;
#        }
#
#        BB[iprev].Rotate_dist(-shift1);
#        BB[i    ].Rotate_dist(-shift2);
#        BB[inext].Rotate_dist(-shift3);
#
#      }
#
#//      Memo1->Lines->Add(AnsiString(i) + ": " + (AnsiString)d12 + "/" + (AnsiString)d23);
#    }
#
#  } // end for it
#
#  Display();
#
#//  Memo1->Clear();
#
#//  Memo1->Lines->Add("Chip:");
#//  Memo1->Lines->Add("Average Pad Distance: " + AnsiString( (int) ((float) W / (float) (N-1))));
#//  Memo1->Lines->Add("Bond pitch: " + AnsiString( (int) PointDistance(BB[0], BB[1])));
#  Memo1->Lines->Add("PCB Bond average: " + AnsiString( (int) d_av ));
#  Memo1->Lines->Add("PCB Bond min    : " + AnsiString( (int) d_min ));
#}

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
if __name__=='__main__':
    bonds = read_chip_pad_definitions('input/Pattern_Top.txt')
    for bond in bonds [:10]: print bond

#void __fastcall TForm1::ButtonLoadClick(TObject *Sender)
#{
#  if (!OpenPatternDialog->Execute()) return;
#  filename = OpenPatternDialog->FileName;
#  std::ifstream f;
#  f.open(filename.c_str(), std::ifstream::in);
#  if (!f.is_open()) {
#    Memo1->Lines->Add("Could not open file!");
#    return;
#  }
#  filename = filename.SubString(0,filename.Length()-4);   // remove '.txt'
#
#  CB_ShiftInterposer->Checked = false;
#
#  // read number of rows and radii
#  f >> Nrow;
#  Memo1->Lines->Add("Nrow = " + AnsiString(Nrow));
#  Ring.resize(Nrow);
#  for (unsigned int i = 0; i<Nrow; i++) {
#    f >> Ring[i];
#    Memo1->Lines->Add("Ring[" + AnsiString(i) + "]: " + AnsiString(Ring[i]));
#  }
#
#  // read angles
#  f >> PHI_L;
#  LE_PhiMin->Text = (AnsiString) PHI_L;
#  PHI_L += 90; PHI_L *= M_PI / 180.0;
#
#  f >> PHI_R;
#  LE_PhiMax->Text = (AnsiString) PHI_R;
#  PHI_R = 90 - PHI_R; PHI_R *= M_PI / 180.0;
#
#  // read pitch
#  f >> pitch;
#  Memo1->Lines->Add("pitch = " + AnsiString(pitch));
#
#  BB.clear();
#  float xpos = 0;
#  do {
#    int readval;
#    f >> readval;
#    if (readval < 0) break;     // stop when we find '-1'
#    Memo1->Lines->Add("Val: " + AnsiString(readval));
#    if (readval > 0) {          // a real bond
#      T_Bond b;
#      b.row  = readval -1;
#      b.l    = Ring[b.row];
#      b.p1.x = xpos;
#      b.p1.y = 0;
#      BB.push_back(b);
#    }
#    xpos += pitch;
#  } while (true);
#  N = BB.size();
#  Memo1->Lines->Add("Number of lines: " + AnsiString(xpos / pitch));
#  Memo1->Lines->Add("Number of bonds: " + AnsiString(N));
#  LE_N->Text = (AnsiString) N;
#
#  float step_phi = (PHI_R - PHI_L) / (N-1);
#  for (unsigned int i = 0; i < N; i++) {
#    float phi  = PHI_L + i * step_phi;
#    BB[i].phi  = phi;
#    BB[i].Calc_p2_from_lphi();
#  }
#
#  Display();
#}

##include <vcl.h>
##include <vector>
##include <cmath>
##include <string>
##include <fstream>
##pragma hdrstop
#
##include "MainUnit.h"
#//---------------------------------------------------------------------------
##pragma package(smart_init)
##pragma resource "*.dfm"
#
#
#// TODO:
#// If several bonds end one inner ring, their spacing could
#// be smaller (e.g. equal to pitch).
#// This could be generalized by adding an information in the input file
#// weather the bond can be minimal or should be spread.
#
#
#TForm1 *Form1;
#
#const float TOGRAD = 180.0 / M_PI;
#
#unsigned int N;                // number of bonds
#unsigned int Nrow;             // number of bond rows
#float        pitch;            // pitch of bond pads
#
#float W, L, PHI_L, PHI_R;
#
#float d_av;                    // average distance of bonds on PCB
#
#struct T_Point2D {
#  float x,y;
#};
#
#
#std::vector<float> Ring;
#AnsiString filename;
#
#std::vector<T_Bond> BB;    // Bonds
#
#float distance (T_Bond & b, T_Point2D & a)
#{
#  float dirx = std::cos(b.phi);
#  float diry = std::sin(b.phi);
#  float lambda = ( (b.p2.x - a.x) * dirx + (b.p2.y-a.y) * diry);
#  float dx = a.x - b.p2.x + lambda * dirx;
#  float dy = a.y - b.p2.y + lambda * diry;
#  return ( lambda > 0.0 ? std::sqrt(dx * dx + dy * dy) : 0.0 );
#}
#
#float PointDistance (T_Bond & B1, T_Bond & B2)  // compare end points
#{
#  return std::max(
#    distance(B1, B2.p2),
#    distance(B2, B1.p2)
#  );
#}
#
#//weg
#__fastcall TForm1::TForm1(TComponent* Owner): TForm(Owner)
#{
#  LE_N->Text = (AnsiString) 49;
#  Nrow = 4;
#  Ring.resize(Nrow);
#  Ring[0] =  700;
#  Ring[1] = 1300;
#  Ring[2] = 1900;
#  Ring[3] = 2500;
#}
#
#//weg
#void TForm1::Display (void)
#{
#  const float BORDER  = 10;
#  float       xmin    = 0;
#  float       xmax    = 0;
#  for (unsigned int i=0; i<N; i++) {
#    xmin = std::min(xmin, BB[i].p2.x);
#    xmax = std::max(xmax, BB[i].p2.x);
#  }
#  float       Wmax    = xmax - xmin;
#  float       Scale   = (Image1->Width - 2 * BORDER) / Wmax;
#  int         YMAX    = Image1->Height - BORDER;
#
#  Image1->Canvas->Brush->Color = clLtGray;
#  Image1->Canvas->FillRect(Rect(0,0,Image1->Width,Image1->Height));
#/*
#  Image1->Canvas->Pen->Color = clRed;
#  Image1->Canvas->Ellipse(DXY-2, YMAX-(DXY-2), DXY+2, YMAX-(DXY+2));
#  Image1->Canvas->Ellipse(DXY+Scale*W-2, YMAX-(DXY-2), DXY+Scale*W+2, YMAX-(DXY+2));
#*/
#  Image1->Canvas->Pen->Color = clBlack;
#  for (unsigned int i=0; i<N; i++) {
#    Image1->Canvas->MoveTo(BORDER+Scale*(BB[i].p1.x-xmin),YMAX-Scale*BB[i].p1.y);
#    Image1->Canvas->LineTo(BORDER+Scale*(BB[i].p2.x-xmin),YMAX-Scale*BB[i].p2.y);
#  }
#}
#
#//weg
#void __fastcall TForm1::ButtonStartClick(TObject *Sender)
#{
#  N      = LE_N         ->Text.ToInt   ();
#  W      = LE_Edge      ->Text.ToDouble();
#  L      = LE_BondLength->Text.ToDouble();
#  PHI_L   = (90.0 + LE_PhiMin->Text.ToDouble()) * M_PI / 180.0;
#  PHI_R   = (90.0 - LE_PhiMax->Text.ToDouble()) * M_PI / 180.0;
#
#  BB.resize(N);
#
#  float step_x   = W / (N-1);
#  float step_phi = (PHI_R-PHI_L) / (N-1);
#  for (unsigned int i = 0; i<N; i++) {
#    BB[i].p1.x = i * step_x;
#    BB[i].p1.y = 0;
#    float phi  = PHI_L + i * step_phi;
#    BB[i].phi  = phi;
#    BB[i].phi  = phi;
#    BB[i].p2.x = BB[i].p1.x + BB[i].l * std::cos(phi);
#    BB[i].p2.y = BB[i].p1.y + BB[i].l * std::sin(phi);
#  }
#
#  Display();
#}
#
#
#//Schritt5
#void __fastcall TForm1::ButtonOutputClick(TObject *Sender)
#{
#  TStringList * l = new TStringList;
#
#  Memo1->Lines->Clear();
#  Memo1->Lines->Add("%!PS-Adobe-3.0 EPSF-3.0");
#  float xmin = 0;
#  float xmax = 0;
#  float ymin = 0;
#  float ymax = 0;
#  for (unsigned int i=0; i<N; i++) {
#    xmin = std::min(xmin, BB[i].p2.x);
#    xmax = std::max(xmax, BB[i].p2.x);
#    ymin = std::min(ymin, BB[i].p2.y);
#    ymax = std::max(ymax, BB[i].p2.y);
#  }
#  const int EDGE = 300;
#
#  int Wbox = xmax - xmin + 2 * EDGE;
#  int Hbox = ymax - ymin + 2 * EDGE;
#  Memo1->Lines->Add("%%BoundingBox: 0 0 " + (AnsiString) Wbox + " " + (AnsiString)Hbox);
#
#  l->LoadFromFile("../FANOUT_Header.ps"); Memo1->Lines->AddStrings(l);
#
#  Memo1->Lines->Add("/DrawSide {");
#  Memo1->Lines->Add("gsave");
#//  Memo1->Lines->Add("-" + (AnsiString) W + " 0 translate");
#  Memo1->Lines->Add("COL_PAD");
#  for (unsigned int i = 0; i < N; i++) {
#    Memo1->Lines->Add(
#      AnsiString((int) ( BB[i].phi / M_PI * 180) ) + " " +
#      AnsiString((int) ( BB[i].p1.x - xmin     ) ) + " " +
#      AnsiString((int)   BB[i].p1.y              ) + " Pad_Small"
#    );
#  }
#  for (unsigned int i = 0; i < N; i++) {
#    Memo1->Lines->Add(
#      AnsiString((int) ( BB[i].phi / M_PI * 180) ) + " " +
#      AnsiString((int) ( BB[i].p2.x - xmin     ) ) + " " +
#      AnsiString((int)   BB[i].p2.y              ) + " Pad_Large"
#    );
#  }
#  Memo1->Lines->Add("COL_BOND");
#  for (unsigned int i = 0; i < N; i++) {
#    Memo1->Lines->Add(
#      AnsiString((int) ( BB[i].p1.x - xmin ) ) + " " +
#      AnsiString((int)   BB[i].p1.y          ) + " moveto " +
#      AnsiString((int) ( BB[i].p2.x - xmin ) ) + " " +
#      AnsiString((int)   BB[i].p2.y          ) + " lineto "
#    );
#  }
#  Memo1->Lines->Add("stroke");
#  Memo1->Lines->Add("grestore");
#  Memo1->Lines->Add("} def");
#
#  Memo1->Lines->Add((AnsiString) EDGE + " dup translate");
#  Memo1->Lines->Add("DrawSide");
#//  Memo1->Lines->Add("1 -1 scale");
#//  Memo1->Lines->Add(" -90 rotate");
#//  Memo1->Lines->Add("DrawSide");
#  Memo1->Lines->Add("showpage");
#
#  Memo1->Lines->SaveToFile(filename + ".ps");
#}
#//---------------------------------------------------------------------------
#
#// Schritt 4
#void __fastcall TForm1::ButtonCreateXYFileClick(TObject *Sender)
#{
#  Memo1->Lines->Clear();
#  Memo1->Lines->Add("#NBonds     : " + AnsiString(N));
#  Memo1->Lines->Add("#Left  Angle: " + AnsiString( (int) (PHI_L * TOGRAD) ));
#  Memo1->Lines->Add("#Right Angle: " + AnsiString( (int) (PHI_R * TOGRAD) ));
#//  Memo1->Lines->Add("#Edge       : " + AnsiString(W));
#  for (unsigned int i=0; i<Nrow; i++)
#    Memo1->Lines->Add("#Ring["+AnsiString(i+1)+"]    : " + AnsiString(Ring[i]));
#  Memo1->Lines->Add("#=======================================#");
#
#  for (unsigned int i=0; i<N; i++) {
#    if ( (BB[i].row == (Nrow-1)) &&  (i%2 == 0) ) {
#      BB[i].l -= 200;
#      BB[i].Calc_p2_from_lphi();
#    }
#    Memo1->Lines->Add(
#    (AnsiString) (int) BB[i].p1.x + " " +
#    (AnsiString) (int) BB[i].p1.y + " " +
#    (AnsiString) (int) (BB[i].phi * TOGRAD - 90) + " " +
#    (AnsiString) (int) BB[i].p2.x + " " +
#    (AnsiString) (int) BB[i].p2.y + " " +
#    (AnsiString) (int) (BB[i].phi * TOGRAD - 90)  );
#  }
#
#  Memo1->Lines->SaveToFile(filename + "_out.txt");
#
#  Display();
#}
#
#
#
#// Schritt 3
#void __fastcall TForm1::ButtonCleanUpClick(TObject *Sender)
#{
#  // find the distances in the OUTER bond row
#  const unsigned int OUTERROW = Nrow-1;
#  std::vector<float> dist;                  // dist[i]: dist between i+1 and i
#  unsigned int i_last = -1;
#
#  float d_bad;
#  unsigned int iproblem = -1;
#  for (unsigned int i = 0; i < N-1; i++) {
#    if (i_last == -1) {                        // no previous bond known
#      if (BB[i].row == OUTERROW) i_last = i;
#    } else {                                   // last bond was in OUTERROW
#      if (BB[i].row == OUTERROW) {             // this as well -> add a distance
#        float dx = BB[i].p2.x - BB[i_last].p2.x;
#        float dy = BB[i].p2.y - BB[i_last].p2.y;
#        d_bad = std::sqrt(dx * dx + dy * dy);
#        // Memo1->Lines->Add(d);
#        if (dx<0) {iproblem = i; break;}
#        i_last = i;
#      }
#    } // if i_last
#  } // for i
#
#  if (iproblem == -1) {
#    Memo1->Lines->Add("Nothing to Clean up");
#  } else {
#    Memo1->Lines->Add("Cleaning up index " + AnsiString(iproblem));
#    Memo1->Lines->Add("d_av = " + AnsiString(d_av));
#    d_bad *= -1;
#    Memo1->Lines->Add("Bad distance = " + AnsiString(d_bad));
#    float shift = (d_av - d_bad) / 2.0;
#    float angle = shift / Ring[OUTERROW];
#
#    for (unsigned int i=0; i<iproblem; i++) {
#      BB[i].Rotate_angle(angle * (float) i / (float) (iproblem-1));
#    }
#
#    for (unsigned int i=iproblem; i<N; i++) {
#      BB[i].Rotate_angle(-angle * (float) (N-1-i) / (float) (N-1-iproblem));
#    }
#
#    Display();
#
#  } // if iproblem
#
#}
#
