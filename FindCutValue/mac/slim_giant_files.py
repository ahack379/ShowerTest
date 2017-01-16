from ROOT import TChain, TFile
import sys

out = TFile("out.root","RECREATE")

print "About to Chain files! "

#get main tree
ch_id = TChain("larlite_id_tree")
ch_mcshr = TChain("mcshower_mcreco_tree")
ch_mctrk = TChain("mctrack_mcreco_tree")
ch_true = TChain("mctruth_generator_tree")
ch_flux = TChain("mcflux_generator_tree")
ch_pot = TChain("potsummary_generator_tree")
ch_gaus = TChain("hit_gaushit_tree")
ch_h02 = TChain("hit_hit02_tree")
ch_vtx = TChain("vertex_numuCC_vertex_tree")
#ch_vtx = TChain("vertex_mcvertex_tree")
ch_trk = TChain("track_pandoraNu_tree")
ch_ass_trk = TChain("ass_pandoraNu_tree")

inputFile = sys.argv[-1]
print "Done chaining files!" 

ch_mcshr.Add(inputFile)
ch_mctrk.Add(inputFile)
ch_true.Add(inputFile)
ch_flux.Add(inputFile) 
ch_gaus.Add(inputFile) 
ch_h02.Add(inputFile) 
ch_vtx.Add(inputFile) 
ch_pot.Add(inputFile) 
ch_id.Add(inputFile) 
ch_trk.Add(inputFile) 
ch_ass_trk.Add(inputFile) 

print "Done adding files!" 

#new tree
ch_new_mcshr = ch_mcshr.CloneTree(0)
ch_new_mctrk = ch_mctrk.CloneTree(0)
ch_new_true = ch_true.CloneTree(0)
ch_new_flux = ch_flux.CloneTree(0)
ch_new_gaus = ch_gaus.CloneTree(0)
ch_new_h02 = ch_h02.CloneTree(0)
ch_new_vtx = ch_vtx.CloneTree(0)
ch_new_pot = ch_pot.CloneTree(0)
ch_new_trk = ch_trk.CloneTree(0)
ch_new_ass_trk = ch_ass_trk.CloneTree(0)
ch_new_id = ch_id.CloneTree(0)

print "Gonna fill mcshower now"

for i in range(ch_mcshr.GetEntries()):
    ch_mcshr.GetEntry(i)
    ch_new_mcshr.Fill()

print "Done filling MCshower"

for i in range(ch_mctrk.GetEntries()):
    ch_mctrk.GetEntry(i)
    ch_new_mctrk.Fill()

print "Done filling MCtrack"

for i in range(ch_true.GetEntries() ):
    ch_true.GetEntry(i)
    ch_new_true.Fill()

print "Done filling MCtrue"

for i in range(ch_flux.GetEntries() ):
    ch_flux.GetEntry(i)
    ch_new_flux.Fill()

for i in range(ch_gaus.GetEntries() ):
    ch_gaus.GetEntry(i)
    ch_new_gaus.Fill()

print "Done filling gaus"

for i in range(ch_h02.GetEntries() ):
    ch_h02.GetEntry(i)
    ch_new_h02.Fill()

for i in range(ch_vtx.GetEntries() ):
    ch_vtx.GetEntry(i)
    ch_new_vtx.Fill()

print "Done filling vertex"

for i in range(ch_pot.GetEntries() ):
    ch_pot.GetEntry(i)
    ch_new_pot.Fill()

for i in range(ch_trk.GetEntries() ):
    ch_trk.GetEntry(i)
    ch_new_trk.Fill()

for i in range(ch_ass_trk.GetEntries() ):
    ch_ass_trk.GetEntry(i)
    ch_new_ass_trk.Fill()

for i in range(ch_id.GetEntries() ):
    ch_id.GetEntry(i)
    ch_new_id.Fill()

print "Done filling!"

ch_new_id.GetCurrentFile().Write()
ch_new_id.GetCurrentFile().Close()

print "DONE!!!!"
