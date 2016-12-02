#ifndef LARLITE_FULLDISTRIBS_CXX
#define LARLITE_FULLDISTRIBS_CXX

#include "FullDistribs.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include "DataFormat/vertex.h"
#include <algorithm>

namespace larlite {

  bool FullDistribs::initialize() {
    
    if( !_energy_tree ){
      _energy_tree = new TTree("energy_tree","");
      _energy_tree->Branch("true_pi0_e",&_true_pi0_e,"true_pi0_e/F");
      _energy_tree->Branch("true_angle",&_true_angle,"true_angle/F");
      _energy_tree->Branch("true_asym",&_true_asym,"true_asym/F");
      _energy_tree->Branch("true_vtx","std::vector<float>",&_true_vtx);
      _energy_tree->Branch("reco_vtx","std::vector<float>",&_reco_vtx);
      _energy_tree->Branch("event",&_event,"event/I");
      }

    if( !_gamma_tree ){
      _gamma_tree = new TTree("gamma_tree","");
      _gamma_tree->Branch("gamma_e",&_gamma_e,"gamma_pi0_e/F");
      _gamma_tree->Branch("rad_l",&_rad_l,"rad_l/F");
      }

    _event = -1;
   std::cout<<"Event list size: "<<_event_list.size() <<std::endl ;

    return true;
  }

  void FullDistribs::Clear(){
    _true_pi0_e = -999;
    _true_angle = -9;
    _true_asym = -9;
    _n_true_pi0 = 0;

    _true_vtx.clear();
    _reco_vtx.clear(); 
    }
  
  bool FullDistribs::analyze(storage_manager* storage) {

    _event ++; 
    Clear() ;

    std::cout<<"\nEvent is : "<<_event <<std::endl ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ){
      std::cout<<"No mctruth..." <<std::endl;
      return false;
      }

    auto parts = ev_mctruth->at(0).GetParticles();
    std::vector<float> start(3,0) ;
    
    for ( auto const & p : parts ){
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){
        _n_true_pi0 ++;
        _true_pi0_e = p.Trajectory().at(0).E()*1000; 
        start[0] = p.Trajectory().at(0).X();
        start[1] = p.Trajectory().at(0).Y();
        start[2] = p.Trajectory().at(0).Z();
        }
      }   

    if ( _n_true_pi0 < 1 ) return false; 

    _true_vtx = start ;

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex"); 
    if(!ev_vtx || !ev_vtx->size() ){
      std::cout<<"No reco vertex..." <<std::endl;
      return false;
      }

    auto & rvtx = ev_vtx->at(0);

    _reco_vtx.emplace_back(rvtx.X());
    _reco_vtx.emplace_back(rvtx.Y());
    _reco_vtx.emplace_back(rvtx.Z());
    
   // Replace pi0 energy with combined shower energy 
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco"); 
    if(!ev_mcs || !ev_mcs->size() ){
      std::cout<<"No mcshower..." <<std::endl;
      return false;
      }

    std::vector<int> shr_ids;
     
    // Find the corresponding mcshowers
    for ( int si = 0; si < ev_mcs->size(); si++){ 

      auto s = ev_mcs->at(si);

      if( s.PdgCode() != 22 || abs(s.MotherPdgCode()) == 13 ) continue; 
      
      auto st = s.Start();
      auto dist = sqrt( pow(st.X() - start[0],2) + pow(st.Y() - start[1],2) + pow(st.Z() - start[2],2) );
      
      if ( dist < 0.001 ){
        shr_ids.emplace_back(si) ;
	_gamma_e = s.Start().E();
	_gamma_tree->Fill();
	}
    }
    
    if( shr_ids.size() == 2 ){

      auto s1 = ev_mcs->at(shr_ids[0]).Start();
      auto s2 = ev_mcs->at(shr_ids[1]).Start();

      auto mag1 = sqrt( s1.Px()*s1.Px()+s1.Py()*s1.Py()+s1.Pz()*s1.Pz() );
      auto mag2 = sqrt( s2.Px()*s2.Px()+s2.Py()*s2.Py()+s2.Pz()*s2.Pz() );
      auto dot = s1.Px()*s2.Px() + s1.Py()*s2.Py() + s1.Pz()*s2.Pz() ;

      float e1 = s1.E() ;
      float e2 = s2.E() ;

      std::cout<<"DOT: "<<dot<<", "<<mag1<<", "<<mag2<<std::endl ; 

      _true_asym = e1 > e2 ? e2/e1 : e1/e2 ;
      _true_angle = acos( dot / mag1 / mag2 ); 
      _energy_tree->Fill();
    }

    return true;
  }

  bool FullDistribs::finalize() {

    if(_fout) { 
      _fout->cd(); 
      _energy_tree->Write(); 
      _gamma_tree->Write(); 
      }
  
    return true;
  }

}
#endif
