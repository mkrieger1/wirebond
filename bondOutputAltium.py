import math

# Delphiscript code from: "Delphiscript Scripts/Pcb/CreateFootprintInLibrary.pas"

#------------------------------------------------------------------------------
# Beginning of Altium .pas file
#------------------------------------------------------------------------------
def altium_header(footprint_name):
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
# Create a pad for each bond
#------------------------------------------------------------------------------
def altium_pad(bond, DX, DY): return """
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
                                  PCBM_BoardRegisteration,NewTrack.I_ObjectAddress);

""" % (bond.p2.x/1000.0, bond.p2.y/1000.0, DX, DY, bond.phi*180/math.pi, bond.name,
       bond.p1.x/1000.0, bond.p1.y/1000.0, bond.p2.x/1000.0, bond.p2.y/1000.0)

#------------------------------------------------------------------------------
# End of Altium .pas file
#------------------------------------------------------------------------------
altium_footer = """
    PCBServer.SendMessageToRobots(CurrentLib.Board.I_ObjectAddress,c_Broadcast,
                                  PCBM_BoardRegisteration,NewPCBLibComp.I_ObjectAddress);
    PCBServer.PostProcess;

    CurrentLib.CurrentComponent := NewPcbLibComp;
    CurrentLib.Board.ViewManager_FullUpdate;
End.
"""


#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
DX = 0.3
DY = 0.15

def bonds_output_altium(bonds, footprint_name, f=None):
    print >> f, altium_header(footprint_name)
    for bond in bonds:
        print >> f, altium_pad(bond, DX, DY)
    print >> f, altium_footer

