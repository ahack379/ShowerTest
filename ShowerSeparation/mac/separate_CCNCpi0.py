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
#my_proc.set_ana_output_file("ana.root");
my_proc.set_output_file("NCpi0.root");

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
aunit = fmwk.SeparateCCNC()
aunit.GetNC(True)

my_proc.add_process(aunit)
my_proc.enable_filter(True)

print
print  "Finished configuring ana_processor. Start event loop!"
print

my_proc.run()


# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
