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
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify analysis output root file name
my_proc.set_ana_output_file("ShowerRecoEff_single_pi0_v2.root");

# Specify data output root file name
#my_proc.set_output_file("qual.root")

clusterer = fmwk.ShowerRecoEff()
my_proc.add_process(clusterer)

#my_proc.set_data_to_write(fmwk.data.kCluster,clusterproducer)
#my_proc.set_data_to_write(fmwk.data.kCluster,'mccluster')
#my_proc.set_data_to_write(fmwk.data.kAssociation,'mccluster')

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()

sys.exit()

