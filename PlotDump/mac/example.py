import sys
from ROOT import gSystem
gSystem.Load("libFindNeutrinos_PlotDump")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing PlotDump..."

sys.exit(0)

