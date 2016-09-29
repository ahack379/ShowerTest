import sys
from ROOT import gSystem
gSystem.Load("libShowerTest_SumCharge")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing SumCharge..."

sys.exit(0)

