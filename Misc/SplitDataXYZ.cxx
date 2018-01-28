#ifndef LARLITE_SPLITDATAXYZ_CXX
#define LARLITE_SPLITDATAXYZ_CXX

#include "SplitDataXYZ.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/potsummary.h"

namespace larlite {

  bool SplitDataXYZ::initialize() {

    if(!_tree){
      _tree = new TTree("tree","tree");
      _tree->Branch("x_pos",&_x_pos,"x_pos/I");
      _tree->Branch("y_pos",&_y_pos,"y_pos/I");
      _tree->Branch("z_pos",&_z_pos,"z_pos/I");
      _tree->Branch("x_neg",&_x_neg,"x_neg/I");
      _tree->Branch("y_neg",&_y_neg,"y_neg/I");
      _tree->Branch("z_neg",&_z_neg,"z_neg/I");
      _tree->Branch("x_pos_POT",&_x_pos_POT,"x_pos_POT/F");
      _tree->Branch("y_pos_POT",&_y_pos_POT,"y_pos_POT/F");
      _tree->Branch("z_pos_POT",&_z_pos_POT,"z_pos_POT/F");
      _tree->Branch("x_neg_POT",&_x_neg_POT,"x_neg_POT/F");
      _tree->Branch("y_neg_POT",&_y_neg_POT,"y_neg_POT/F");
      _tree->Branch("z_neg_POT",&_z_neg_POT,"z_neg_POT/F");
    }

    _x_pos = 0;
    _y_pos = 0;
    _z_pos = 0;
    _x_neg = 0;
    _y_neg = 0;
    _z_neg = 0;

    _x_pos_POT = 0;
    _y_pos_POT = 0;
    _z_pos_POT = 0;
    _x_neg_POT = 0;
    _y_neg_POT = 0;
    _z_neg_POT = 0;

    return true;
  }
  
  bool SplitDataXYZ::analyze(storage_manager* storage) {
  

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
    }

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();

    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    // Now calculate the total POT + total numu neutrinos 
    auto ev_pot = storage->get_subrundata<potsummary>("generator"); 

    if ( xyz[0] > 256.35 / 2 ) {
      _x_pos ++ ;
      if( storage->subrun_id() != storage->last_subrun_id() )
        _x_pos_POT += ev_pot->totgoodpot ;
    }
    else {
       _x_neg ++;
      if( storage->subrun_id() != storage->last_subrun_id() )
        _x_neg_POT += ev_pot->totgoodpot ;
    }
     
    if ( xyz[1] > 233. / 2 ){
      _y_pos ++;
      if( storage->subrun_id() != storage->last_subrun_id() )
        _y_pos_POT += ev_pot->totgoodpot ;
    }
    else{
      _y_neg ++;
      if( storage->subrun_id() != storage->last_subrun_id() )
        _y_neg_POT += ev_pot->totgoodpot ;
    }

    if ( xyz[2] > 1036.8 / 2 ){
       _z_pos ++;
      if( storage->subrun_id() != storage->last_subrun_id() )
        _z_pos_POT += ev_pot->totgoodpot ;
    }
    else{
       _z_neg ++;
      if( storage->subrun_id() != storage->last_subrun_id() )
        _z_neg_POT += ev_pot->totgoodpot ;
    }

    _tree->Fill();

    return true;
  }

  bool SplitDataXYZ::finalize() {
    
    if(_fout) { _fout->cd(); _tree->Write(); }
  
    return true;
  }

}
#endif
