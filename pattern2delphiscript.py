#!/usr/bin/python

# Delphiscript code from: "Delphiscript Scripts/Pcb/CreateFootprintInLibrary.pas"

#------------------------------------------------------------------------------
# begin of .pas file
#------------------------------------------------------------------------------
def header(footprint_name):
    return """Var
    CurrentLib    : IPCB_Library;
    NewPCBLibComp : IPCB_LibComponent;
    NewPad        : IPCB_Pad;
    NewTrack      : IPCB_Track;

Begin
    If PCBServer = Nil Then Exit;
    CurrentLib := PcbServer.GetCurrentPCBLibrary;
    If CurrentLib = Nil Then Exit;

    NewPCBLibComp := PCBServer.CreatePCBLibComp;
    NewPcbLibComp.Name := '%s';

    CurrentLib.RegisterComponent(NewPCBLibComp);

    PCBServer.PreProcess;
""" % footprint_name


#------------------------------------------------------------------------------
# repeat for every pad
#------------------------------------------------------------------------------
def bond(x1, y1, x2, y2, dx, dy, angle, name): return """
    NewPad := PcbServer.PCBObjectFactory(ePadObject,eNoDimension,eCreate_Default);
    NewPad.X        := MMsToCoord(%.3f);
    NewPad.Y        := MMSToCoord(%.3f);
    NewPad.TopXSize := MMsToCoord(%.3f);
    NewPad.TopYSize := MMsToCoord(%.3f);   
    NewPad.HoleSize := MMsToCoord(0);
    NewPad.Rotation := %i;
    NewPad.Layer    := eTopLayer;
    NewPad.Name     := '%s';
    NewPad.TopShape := eRectangular;
    NewPCBLibComp.AddPCBObject(NewPad);
    PCBServer.SendMessageToRobots(NewPCBLibComp.I_ObjectAddress,c_Broadcast,
                                  PCBM_BoardRegisteration,NewPad.I_ObjectAddress);
    NewTrack := PcbServer.PCBObjectFactory(eTrackObject,eNoDimension,eCreate_Default);
    NewTrack.X1 := MMsToCoord(%.3f);
    NewTrack.Y1 := MMsToCoord(%.3f);
    NewTrack.X2 := MMsToCoord(%.3f);
    NewTrack.Y2 := MMsToCoord(%.3f);
    NewTrack.Layer := eMechanical15;
    NewTrack.Width := MMsToCoord(0.01);
    NewPCBLibComp.AddPCBObject(NewTrack);
    PCBServer.SendMessageToRobots(NewPCBLibComp.I_ObjectAddress,c_Broadcast,
                                  PCBM_BoardRegisteration,NewTrack.I_ObjectAddress)

""" % (x2, y2, dx, dy, angle, name, x1, y1, x2, y2)

#------------------------------------------------------------------------------
# end of .pas file
#------------------------------------------------------------------------------
footer = """
    PCBServer.SendMessageToRobots(CurrentLib.Board.I_ObjectAddress,c_Broadcast,
                                  PCBM_BoardRegisteration,NewPCBLibComp.I_ObjectAddress);
    PCBServer.PostProcess;

    CurrentLib.CurrentComponent := NewPcbLibComp;
    CurrentLib.Board.ViewManager_FullUpdate;
End.
"""


#------------------------------------------------------------------------------
# read & process pattern file
#------------------------------------------------------------------------------
def read_pattern_file(pattern_file, angle_offset=0, pad_number_offset=0):
    s = ''
    with open(pattern_file) as f:
        for (i, line) in enumerate(f):
            if not line.startswith('#'):
                l = line.split()
                [x1, y1] = [float(x)/1000 for x in l[0:2]]
                [x2, y2] = [float(x)/1000 for x in l[3:5]]
                angle = int(l[5])

                s += bond(x1, y1, x2, y2, DX, DY,
                          angle+angle_offset, str(i+pad_number_offset))
    return s


#------------------------------------------------------------------------------
# pad geometry definition
#------------------------------------------------------------------------------
angle_offset = {'E': 0, 'N': 90, 'W': 180, 'S': 270}
DX = 0.3
DY = 0.15


#------------------------------------------------------------------------------
# main
#------------------------------------------------------------------------------
def create_delphiscript_file(footprint_name, pattern_file):
    print header(footprint_name)
    #print read_pattern_file(pattern_file, angle_offset['E']) # "north" edge of the chip
    print read_pattern_file(pattern_file, angle_offset['N']) # "north" edge of the chip
    #print read_pattern_file(pattern_file, angle_offset['W']) # "north" edge of the chip
    #print read_pattern_file(pattern_file, angle_offset['S']) # "north" edge of the chip
    print footer

if __name__=='__main__':
    import sys

    if len(sys.argv) < 3:
        print 'usage: %s <footprint_name> <pattern_filename>' % sys.argv[0]

    else:
        footprint_name = sys.argv[1]
        pattern_file = sys.argv[2]
        create_delphiscript_file(footprint_name, pattern_file)

