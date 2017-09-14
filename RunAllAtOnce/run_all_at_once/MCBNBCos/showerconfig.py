import sys

# Import the needed modules.  Need larlite and several of it's namespaces
from ROOT import gSystem,TMath
from larlite import larlite as fmwk
from larlite import larutil
from recotool import cmtool, showerreco

def getShowerRecoAlgModular():
    
    alg = showerreco.ShowerRecoAlgModular()
    alg.SetDebug(False)
    alg.SetVerbose(False)
    
    angle3D = showerreco.Angle3DFormula()
    angle3D.setMaxAngleError(0.1)
    angle3D.setValidateDirection(True)
    angle3D.setVerbosity(False)

    angle3D = showerreco.Angle3DFromVtxQweighted()
    
    energy = showerreco.LinearEnergy()
    energy.SetElectronLifetime(1e20) 
    energy.SetRecombFactor(0.577) 
    #energy.SetElecGain(243) # e- / ADC DATA value, MCC8 Data
    energy.SetElecGain(198) # MCC8 MC value
    energy.setVerbosity(False)
    energy.SetFillTree(True)

    dqdx = showerreco.dQdxModule()
    
    dedx = showerreco.dEdxFromdQdx()
    dedx.SetUsePitch(False)
    dedx.setVerbosity(True)

    geoModule = showerreco.GeoModule()
    geoModule.setFlipShowerDirection(False)
    
    alg.AddShowerRecoModule( showerreco.StartPoint3DModule() )
    alg.AddShowerRecoModule(angle3D)
    alg.AddShowerRecoModule(energy)
    alg.AddShowerRecoModule( geoModule )
    
    alg.PrintModuleList()
    
    return alg


