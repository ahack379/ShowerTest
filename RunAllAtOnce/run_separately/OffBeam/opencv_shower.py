from showerconfig import getShowerRecoAlgModular 

import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk
from ROOT import protoshower

def DefaultShowerReco3D():

    ana_unit = fmwk.ShowerReco3D()
    ana_unit.SetRequirePDG11(False)
    sralg = getShowerRecoAlgModular()
    ana_unit.AddShowerAlgo(sralg)

    return ana_unit

my_proc = fmwk.ana_processor()

for x in xrange(len(sys.argv)-1 ): 
    my_proc.add_input_file(sys.argv[x+1])

cfg  = 'MCBNB_cosmics_BBox.fcl'
#name = sys.argv[1]

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

my_proc.enable_filter(False)

# Specify data output root file name
my_proc.set_output_file("opencv_shower_out.root") #%s_output.root" % name); 

my_proc.set_ana_output_file("opencv_ana.root");

myunit = fmwk.LArImageHit()
myunit.set_config(cfg)
my_proc.add_process(myunit,False)

protoshoweralg = protoshower.ProtoShowerAlgOpenCV()
shr_unit = DefaultShowerReco3D()
shr_unit.GetProtoShowerHelper().setProtoShowerAlg( protoshoweralg )
shr_unit.SetInputProducer("ImageClusterHit")
shr_unit.SetOutputProducer("showerreco")
my_proc.add_process(shr_unit)

# Output we save
# HitRemoval Stuff
my_proc.set_data_to_write(fmwk.data.kHit,         "gaushit"        )
my_proc.set_data_to_write(fmwk.data.kCluster,     "pandoraCosmic"  )
my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraCosmic" )

# OpenCV Stuff
my_proc.set_data_to_write(fmwk.data.kCluster,     "ImageClusterHit"  )
my_proc.set_data_to_write(fmwk.data.kAssociation, "ImageClusterHit"  )
my_proc.set_data_to_write(fmwk.data.kPFParticle,  "ImageClusterHit"  )

# Showerreco Stuff
my_proc.set_data_to_write(fmwk.data.kShower,      "showerreco")
my_proc.set_data_to_write(fmwk.data.kAssociation, "showerreco")

# Other Things 
my_proc.set_data_to_write(fmwk.data.kVertex,  "numuCC_vertex"  )
my_proc.set_data_to_write(fmwk.data.kTrack,   "numuCC_track"  )
my_proc.set_data_to_write(fmwk.data.kTrack,   "pandoraNu"  )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run() 

sys.exit()

