from algoconfig import loadAlgo
from algoconfig_nu import loadAlgo as loadAlgo_nu
from showerconfig import getShowerRecoAlgModular 

import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from larlite import larlite as fmwk

#from ROOT import twodimtools
from ROOT import protoshower

#a = twodimtools.Linearity()
#print a

def DefaultShowerReco3D():

    ana_unit = fmwk.ShowerReco3D()
    ana_unit.SetRequirePDG11(False)
    sralg = getShowerRecoAlgModular()
    ana_unit.AddShowerAlgo(sralg)

    return ana_unit

my_proc = fmwk.ana_processor()

for x in xrange(len(sys.argv)-1 - 1): 
    print sys.argv[x+2]
    my_proc.add_input_file(sys.argv[x+2])

cfg  = 'MCBNB_cosmics_BBox.fcl'
name = sys.argv[1]

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

my_proc.enable_filter(True)

# Specify data output root file name
my_proc.set_output_file("%s_output.root" % name); 

my_proc.set_ana_output_file("analysis.root");

#my_proc.enable_event_alignment(False)

# prepare the various hit removal stages
my_proc.add_process( loadAlgo("ROIRemoval") )
my_proc.add_process( loadAlgo("VertexSlopeCorrelation") )
my_proc.add_process( loadAlgo("VertexAngleCorrelation") )
my_proc.add_process( loadAlgo("RemoveDeltaRays") )

my_proc.add_process(loadAlgo_nu("VertexProximityRemoval") )
my_proc.add_process(loadAlgo_nu("PandoraLinearRemoval") )
my_proc.add_process(loadAlgo_nu("TrackDeltaRayRemoval") )
my_proc.add_process(loadAlgo_nu("SimpleClusterer") )
my_proc.add_process(loadAlgo_nu("VertexTrackRemoval") )
my_proc.add_process(loadAlgo_nu("LinearRemoval") )
my_proc.add_process(loadAlgo_nu("RemoveHitsNearVtx") )

ratio = fmwk.DistanceCut()
ratio.SetRatioCut(0.22) 
ratio.SetRadius(60)
my_proc.add_process(ratio,True)

myunit = fmwk.LArImageHit()
myunit.set_config(cfg)
my_proc.add_process(myunit,True)

protoshoweralg = protoshower.ProtoShowerAlgOpenCV()
shr_unit = DefaultShowerReco3D()
shr_unit.GetProtoShowerHelper().setProtoShowerAlg( protoshoweralg )
shr_unit.SetInputProducer("ImageClusterHit")
shr_unit.SetOutputProducer("showerreco")
my_proc.add_process(shr_unit)

pi0 = fmwk.Pi0Cuts()
my_proc.add_process(pi0,True)

# Output we save
# HitRemoval Stuff
my_proc.set_data_to_write(fmwk.data.kHit,         "gaushit"        )
my_proc.set_data_to_write(fmwk.data.kCluster,     "pandoraCosmic"  )
my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraCosmic" )
#my_proc.set_data_to_write(fmwk.data.kCluster,     "sc"            )   
#my_proc.set_data_to_write(fmwk.data.kAssociation, "sc"            )   

# OpenCV Stuff
my_proc.set_data_to_write(fmwk.data.kCluster,     "ImageClusterHit"  )
my_proc.set_data_to_write(fmwk.data.kAssociation, "ImageClusterHit"  )
my_proc.set_data_to_write(fmwk.data.kPFParticle,  "ImageClusterHit"  )

# Showerreco Stuff
my_proc.set_data_to_write(fmwk.data.kShower,      "showerreco")
my_proc.set_data_to_write(fmwk.data.kAssociation, "showerreco")

# Pi0Reco Stuff
my_proc.set_data_to_write(fmwk.data.kShower,      "pi0_candidate_showers")
my_proc.set_data_to_write(fmwk.data.kAssociation, "pi0_candidate_showers")

# Other Things 
my_proc.set_data_to_write(fmwk.data.kVertex,  "numuCC_vertex"  )
my_proc.set_data_to_write(fmwk.data.kTrack,   "numuCC_track"  )
my_proc.set_data_to_write(fmwk.data.kTrack,   "pandoraNu"  )

# MC Stuff 
my_proc.set_data_to_write(fmwk.data.kMCShower,    "mcreco"         )
my_proc.set_data_to_write(fmwk.data.kMCTrack,     "mcreco"         )
my_proc.set_data_to_write(fmwk.data.kMCTruth,     "generator"      )
#my_proc.set_data_to_write(fmwk.data.kCluster,  "mccluster"  )
#my_proc.set_data_to_write(fmwk.data.kAssociation,  "mccluster"  )
my_proc.set_data_to_write(fmwk.data.kVertex,      "mcvertex"       )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run(0,20) 

sys.exit()

