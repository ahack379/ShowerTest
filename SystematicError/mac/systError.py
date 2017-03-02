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
my_proc.set_ana_output_file("from_test_ana_you_can_remove_me.root");

ana = fmwk.RepoSelectionII()
# SelII Variables
ana.SetBounds_DistToEdgeXYZ(5.,5.,0., 38.,38.,20.)
ana.SetBounds_PEThresh           (40., 60.)
ana.SetBounds_Trk2FlashDist      (55., 85)
ana.SetBounds_MinTrk2VtxDist     (0., 7.)
ana.SetBounds_MinTrackLen        (5., 50.)
ana.SetBounds_MaxCosineAngle     (0.85, 1.0)
ana.SetBounds_MaxCosy1stTrk      (0.5, 0.85)
ana.SetBounds_MinTrackLen2ndTrk  (20., 50.)
ana.SetBounds_MaxCosySingle      (0.5, 0.7)
ana.SetBounds_MinTrackLenSingle  (30., 40.)
ana.SetBounds_MindEdxRatioSingle (1.0, 2. )
ana.SetBounds_MaxTrkLengthySingle(15., 45.)
ana.SetBounds_MinStartdEdx1stTrk (2.1, 3.5)
ana.SetBounds_MaxEnddEdx1stTrk   (1.0, 10.)

# HitRatio Variables
ana.SetBounds_Radius        (40., 75.)
ana.SetBounds_HitRatio      (0.15,0.38)
ana.SetBounds_MinHitsRequired (10,50)

# ShowerPair Variables
ana.SetBounds_OpeningAngle  (0.3, 0.65)
ana.SetBounds_MaxIP         (2., 8.)
ana.SetBounds_MaxRadLength  (40, 80)

# Other Variabels 
ana.SetIters(50)

my_proc.add_process(ana)

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run(0,10);

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
