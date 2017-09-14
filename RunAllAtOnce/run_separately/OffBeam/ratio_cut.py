import sys
from larlite import larlite as fmwk

my_proc = fmwk.ana_processor()

for x in xrange(len(sys.argv)-1): 
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

my_proc.enable_filter(True)

# Specify data output root file name
my_proc.set_output_file("ratiocut_out.root"); 

my_proc.set_ana_output_file("ratiocut_ana.root");

ratio = fmwk.DistanceCut()
ratio.SetRatioCut(0.19)#0.22) 
ratio.SetRadius(50)#60)
ratio.SetVtxProducer("numuCC_vertex")
my_proc.add_process(ratio,True)

# Output we save
# HitRemoval Stuff
my_proc.set_data_to_write(fmwk.data.kHit,         "gaushit"        )
my_proc.set_data_to_write(fmwk.data.kCluster,     "pandoraCosmic"  )
my_proc.set_data_to_write(fmwk.data.kAssociation, "pandoraCosmic" )
#my_proc.set_data_to_write(fmwk.data.kCluster,     "sc"            )   
#my_proc.set_data_to_write(fmwk.data.kAssociation, "sc"            )   

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
#my_proc.set_data_to_write(fmwk.data.kVertex,      "mcvertex"       )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run() 

sys.exit()

