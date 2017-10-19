#ifndef LARLITE_NEUTRINOMODE_CXX
#define LARLITE_NEUTRINOMODE_CXX

#include "NeutrinoMode.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/track.h"
#include "DataFormat/shower.h"
#include "DataFormat/vertex.h"

#include "LArUtil/GeometryHelper.h"

namespace larlite {

  bool NeutrinoMode::initialize() {    

    _event = -1; 

    if ( !_tree ){
      _tree = new TTree("tree","");
      _tree->Branch("event",&_event,"event/I");
      _tree->Branch("bkgd_id",&_bkgd_id,"bkgd_id/I");
      _tree->Branch("nu_mode",&_nu_mode,"nu_mode/I");

   }

    return true;
  }

  void NeutrinoMode::clear(){

    _bkgd_id = -1 ;
    _nu_mode = -1 ;

  }
  
  bool NeutrinoMode::analyze(storage_manager* storage) {

    _event++ ;
    clear();

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
    }

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();
    _nu_mode = nu.Mode();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();
    auto e = traj.at(traj.size() - 1).E();

    bool infv = true;
    if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
      infv = false;

    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0;

    for ( auto const & p : parts ){
      if( p.StatusCode() == 1 && p.PdgCode() == 111 )
        n_pi0 ++;
    }   
    
    if( nu.Nu().PdgCode() == 14 && nu.CCNC() == 0 && n_pi0 == 1 && infv && e > 0.275 ) 
      _bkgd_id = 2;
    else if( nu.CCNC() == 0 && n_pi0 == 0 ) 
      _bkgd_id = 3;
    else if( nu.CCNC() == 1 && n_pi0 > 0 ) 
      _bkgd_id = 4;
    else if( nu.CCNC() == 1 && n_pi0 == 0 ) 
      _bkgd_id = 5;
    else 
      _bkgd_id = 6;
    
    _tree->Fill();    

    return true;
  }

  bool NeutrinoMode::finalize() {

    if ( _fout ){
      _fout->cd();
      _tree->Write();
    }
  
    return true;
  }

}
#endif
