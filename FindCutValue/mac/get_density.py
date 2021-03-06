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
my_proc.set_io_mode(fmwk.storage_manager.kREAD)


# Specify output root file name
my_proc.set_ana_output_file("density_ana.root") #sys.argv[-1])#"ana.root");

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
ana = fmwk.VtxDensity()
ana.GetPi0s(False)
my_proc.add_process(ana)

####
#my_proc.set_output_file("DENSITY_TEST.root") #sys.argv[-1])#"ana.root");

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run()#0,5050);

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
