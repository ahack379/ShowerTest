#ifndef LARLITE_RLENERGYCUTSTUDY_CXX
#define LARLITE_RLENERGYCUTSTUDY_CXX

#include "RLEnergyCutStudy.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"

namespace larlite {

  bool RLEnergyCutStudy::initialize() {

    if( !_shower_tree ){
     _shower_tree = new TTree("shower_tree","");
     _shower_tree->Branch("_event",&_event,"event/I");
     _shower_tree->Branch("_signal",&_signal,"signal/I");
     _shower_tree->Branch("_shower_e",&_shower_e,"shower_e/F");
     _shower_tree->Branch("_shower_rl",&_shower_rl,"shower_rl/F");
    }

    _event = -1;
    _signal = false; 

    return true;
  }

 void RLEnergyCutStudy::clear(){

    _shower_rl = -10;
    _shower_e  = -10;
    _signal = false ;

  }

  bool RLEnergyCutStudy::analyze(storage_manager* storage) {

    _event++; 
    clear();

    std::cout<<"\nEvent : "<<_event <<std::endl;

    auto ev_mctruth = storage->get_data<event_mctruth>("generator");

    if( !ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Not enough mctruth..." <<ev_mctruth->size()<<std::endl; return false; }

    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");

    if( !ev_mcs || !ev_mcs->size() ) {
      std::cout<<"Not enough mcshowers..." <<ev_mcs->size()<<std::endl; return false; }

    if( !ev_mcs || !ev_mcs->size() ) {
      std::cout<<"Not enough mcshowers..." <<ev_mcs->size()<<std::endl; return false; }


    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();
    auto e = traj.at(traj.size() - 1).E();

    if( xyz[0] < 0 || xyz[0] > 256.35 || xyz[1] > 116.5 || xyz[1] < -116.5 || xyz[2] < 0 || xyz[2] > 1036.8 )
            return false;

    int max_it = -1;
    float max_e = 0.; 

    auto parts = truth.GetParticles() ;

    int n_pi0 = 0;
    int n_mu = 0;

    for ( auto const & p : parts ){
      if ( p.PdgCode() == 111 && p.StatusCode() == 1 )
        n_pi0++;

      if ( p.PdgCode() == 13 && p.StatusCode() == 1 )
        n_mu++;
    }

    if ( n_pi0 == 1 && n_mu == 1 ) _signal = true;

    for ( int i = 0; i < ev_mcs->size(); i++ ){ 

      auto s = ev_mcs->at(i); 
      auto st = s.Start() ;

      auto dist = sqrt( pow(st.X() - xyz[0],2) + pow(st.Y() - xyz[1],2) + pow(st.Z() - xyz[2],2) ); 
    
      if ( dist < 0.001 && s.DetProfile().E() > max_e ) {
        std::cout<<"energy : "<<s.DetProfile().E() <<std::endl ;

        auto det = s.DetProfile() ;

        max_e = det.E() ;
        max_it = i ;

        auto dist_det = sqrt( pow(det.X() - xyz[0],2) + pow(det.Y() - xyz[1],2) + pow(det.Z() - xyz[2],2) ); 

        _shower_e = max_e;
        _shower_rl = dist_det ;

	//std::cout<<"dist_det : "<<dist_det <<", "<<det.X()<<", "<<det.Y()<<", "<<det.Z()<<std::endl ; 

      }
    }  

    _shower_tree->Fill();


    return true;
  }

  bool RLEnergyCutStudy::finalize() {

    if(_fout) {
      _fout->cd();
      _shower_tree->Write();
    }
  
    return true;
  }

}
#endif
