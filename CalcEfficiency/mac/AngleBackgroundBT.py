import sys
if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)


from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    fname = sys.argv[x+1]
    my_proc.add_input_file(fname)
    
# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

search = fmwk.SaveSel2VtxTrk()
search.SetVerbose(False)
my_proc.add_process(search)

clusterer = fmwk.MCNuClusterer()
clusterer.setMinEnergy(5) #MeV
my_proc.add_process(clusterer)

# Specify data output root file name
my_proc.set_ana_output_file("sel2_noMIPangle_ana.root")
my_proc.set_output_file("/Volumes/UBooNEDataMCC8/MCC8_3/OutputFiles/MCBNBCos/anafiles/sel2_noMIPangle_out.root")

clusterer = fmwk.AngleBackgroundBT()
clusterer.UseMCSample(True)
clusterer.GetPi0Info(False)
clusterer.GetSingleShowerInfo(False)
my_proc.add_process(clusterer)

# Backtracker
my_proc.set_data_to_write(fmwk.data.kCluster,     "mccluster")
my_proc.set_data_to_write(fmwk.data.kAssociation, "mccluster")

# Other Things 
my_proc.set_data_to_write(fmwk.data.kVertex,  "numuCC_vertex"  )
my_proc.set_data_to_write(fmwk.data.kTrack,   "numuCC_track"  )
my_proc.set_data_to_write(fmwk.data.kVertex,   "mcvertex"  )

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

