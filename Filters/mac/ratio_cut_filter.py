import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE(s)\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify output root file name
#my_proc.set_ana_output_file("density_ana.root") #sys.argv[-1])#"ana.root");
#my_proc.set_output_file("/Volumes/UBooNEData/mcc8/cosmic_bnb/anafiles/mcc8_bnbcos_post_ratiocut.root")
my_proc.set_output_file("/Volumes/UBooNEData/mcc8/cosmic_bnb/anafiles/mcc8_bnbcos_post_ratiocut_wpandora.root")

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
ratio = fmwk.RatioCut()
ratio.SetRatioCut(0.225) #1875) #0.195) #24)
my_proc.add_process(ratio)

my_proc.enable_filter(True)

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.set_data_to_write(fmwk.data.kMCTruth,'generator')
my_proc.set_data_to_write(fmwk.data.kMCTrack,'mcreco')
my_proc.set_data_to_write(fmwk.data.kMCShower,'mcreco')
my_proc.set_data_to_write(fmwk.data.kVertex,'numuCC_vertex')
my_proc.set_data_to_write(fmwk.data.kTrack,'pandoraNu')
my_proc.set_data_to_write(fmwk.data.kTrack,'numuCC_track')
my_proc.set_data_to_write(fmwk.data.kHit,'gaushit')
#my_proc.set_data_to_write(fmwk.data.kPOTSummmary,'generator')

# Let's run it.
my_proc.run() #5050)

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
