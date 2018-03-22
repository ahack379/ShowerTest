#ifndef LARLITE_MCVARIABLESTUDY_CXX
#define LARLITE_MCVARIABLESTUDY_CXX

#include "MCVariableStudy.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/shower.h"
#include <algorithm>

namespace larlite {

  bool MCVariableStudy::initialize() {
    
    if( !_pi0_tree ){
      _pi0_tree = new TTree("pi0_tree","");
      _pi0_tree->Branch("run_id",&_run_id,"run_id/F");
      _pi0_tree->Branch("subrun_id",&_subrun_id,"subrun_id/F");
      _pi0_tree->Branch("event_id",&_event_id,"event_id/F");
      _pi0_tree->Branch("true_pi0_e",&_true_pi0_e,"true_pi0_e/F");
      _pi0_tree->Branch("true_angle",&_true_angle,"true_angle/F");
      _pi0_tree->Branch("true_asym",&_true_asym,"true_asym/F");
      _pi0_tree->Branch("reco_pi0_e",&_reco_pi0_e,"reco_pi0_e/F");
      _pi0_tree->Branch("true_pi0_mom",&_true_pi0_mom,"true_pi0_mom/F");
      _pi0_tree->Branch("true_gamma_e_min",&_true_gamma_e_min,"true_gamma_e_min/F");
      _pi0_tree->Branch("true_gamma_e_max",&_true_gamma_e_max,"true_gamma_e_max/F");
      _pi0_tree->Branch("true_RL_maxE",&_true_RL_maxE,"true_RL_maxE/F");
      _pi0_tree->Branch("true_RL_minE",&_true_RL_minE,"true_RL_minE/F");
      _pi0_tree->Branch("true_nu_e",&_true_nu_e,"true_nu_e/F");
      _pi0_tree->Branch("true_mu_mom",&_true_mu_mom,"true_mu_mom/F");
      _pi0_tree->Branch("true_mu_len",&_true_mu_len,"true_mu_len/F");
      _pi0_tree->Branch("true_mu_theta",&_true_mu_theta,"true_mu_theta/F");
      _pi0_tree->Branch("true_mu_phi",&_true_mu_phi,"true_mu_phi/F");
      _pi0_tree->Branch("event",&_event,"event/I");
      }

    if( !_gamma_tree ){
      _gamma_tree = new TTree("gamma_tree","");
      _gamma_tree->Branch("true_gamma_e",&_true_gamma_e,"true_gamma_pi0_e/F");
      _gamma_tree->Branch("reco_gamma_e",&_reco_gamma_e,"reco_gamma_pi0_e/F");
      _gamma_tree->Branch("true_reco_dot",&_true_reco_dot,"true_reco_dot/F");
      _gamma_tree->Branch("true_rad_l",&_true_rad_l,"true_rad_l/F");
      _gamma_tree->Branch("reco_rad_l",&_reco_rad_l,"reco_rad_l/F");
      }

    _event = -1;

   std::cout<<"Event list size: "<<_event_list.size() <<std::endl ;

   _n_signals = 0;

    return true;
  }

  void MCVariableStudy::Clear(){
    _true_pi0_e = -999;
    _true_angle = -9;
    _true_asym = -9;
    _reco_pi0_e = -999;
    _true_gamma_e_max = -999 ;
    _true_gamma_e_min = -999 ;
    _true_pi0_mom = -999;
    _true_nu_e = -999;
    _true_RL_maxE = -999;
    _true_RL_minE = -999;

    _true_mu_mom = -999;
    _true_mu_len = -999;
    _true_mu_theta = -999;
    _true_mu_phi = -999;

    _n_true_pi0 = 0;

    _true_gamma_e = -9;
    _reco_gamma_e = -9;
    _true_rad_l = -9;
    _reco_rad_l = -9;
    _true_reco_dot = -9;
    }
  
  bool MCVariableStudy::analyze(storage_manager* storage) {

    _event ++; 
    Clear() ;

    _run_id = storage->run_id();
    _subrun_id = storage->subrun_id();
    _event_id = storage->event_id();

    //std::cout<<"\nEvent is : "<<_event <<std::endl ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ){
      std::cout<<"No mctruth..." <<std::endl;
      return false;
    }

    auto & nu  = ev_mctruth->at(0).GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    bool infv = true;

    if( xyz[0] < 20 || xyz[0] > 236.25 || xyz[1] < -96.5 || xyz[1] > 96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
      infv = false ;

    auto parts = ev_mctruth->at(0).GetParticles();
    int n_mu = 0;
    int n_pi0 = 0;
    
    std::vector<float> start(3,0) ;
    
    for ( auto const & p : parts ){
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){
        _n_true_pi0 ++;
	n_pi0++;
        _true_pi0_e = p.Trajectory().at(0).E()*1000; 
        start[0] = p.Trajectory().at(0).X();
        start[1] = p.Trajectory().at(0).Y();
        start[2] = p.Trajectory().at(0).Z();
	_true_pi0_mom = sqrt( pow(p.Trajectory().at(0).Px(),2) + pow(p.Trajectory().at(0).Py(),2) + pow(p.Trajectory().at(0).Pz(),2) )*1000; 
      }

      if( p.StatusCode() == 1 && p.PdgCode() == 13 )
        n_mu ++; 
    }

    //if( n_mu != 1 || _n_true_pi0 != 1) return false ; 

    if( nu.Nu().PdgCode() != 14 || n_pi0 != 1 || nu.CCNC() != 0 || !infv ) return false;

    _n_signals++ ;


    _true_nu_e = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().at(0).E() ;
      
   // Replace pi0 energy with combined shower energy 
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco"); 
    if(!ev_mcs || !ev_mcs->size() ){
      std::cout<<"No mcshower..." <<std::endl;
      return false;
    }

    std::vector<int> shr_ids;
    std::map<float,int> shr_map ;
     
    for ( int si = 0; si < ev_mcs->size(); si++){ 

      auto s = ev_mcs->at(si);

      if( s.PdgCode() != 22 || s.MotherPdgCode() != 111 || s.Origin() != 1 ) continue; //|| abs(s.MotherPdgCode()) == 13 ) continue;
      
      auto st = s.Start();
      auto dist = sqrt( pow(st.X() - start[0],2) + pow(st.Y() - start[1],2) + pow(st.Z() - start[2],2) );
      
      if ( dist < 0.0001 ){
        shr_ids.emplace_back(si) ;
        shr_map.emplace(1./st.E(), si);
        //std::cout<<"Dist: "<<dist<<", energies shr: "<<ev_mcs->at(si).Start().E()<<std::endl;
	  }
    }
    
    if( shr_ids.size() > 1 ){ // == 2){

      int leading_shower_id = -1;
      int subleading_shower_id = -1;

      int ii = 0;
      for ( auto const & map_element : shr_map ){
        if( ii == 0 ){
          leading_shower_id = map_element.second;
          ii++;
        }
        if( ii == 1 ){
          subleading_shower_id = map_element.second;
          break;
        }
      }

      auto s1 = ev_mcs->at(leading_shower_id).Start();
      auto s2 = ev_mcs->at(subleading_shower_id).Start();

      auto d1 = ev_mcs->at(leading_shower_id).DetProfile();
      auto d2 = ev_mcs->at(subleading_shower_id).DetProfile();

      auto mag1 = sqrt( s1.Px()*s1.Px()+s1.Py()*s1.Py()+s1.Pz()*s1.Pz() );
      auto mag2 = sqrt( s2.Px()*s2.Px()+s2.Py()*s2.Py()+s2.Pz()*s2.Pz() );
      auto dot = s1.Px()*s2.Px() + s1.Py()*s2.Py() + s1.Pz()*s2.Pz() ;

      auto radL1 = sqrt( pow(d1.X() - xyz[0],2) + pow(d1.Y() - xyz[1],2) + pow(d1.Z() - xyz[2],2) ); 
      auto radL2 = sqrt( pow(d2.X() - xyz[0],2) + pow(d2.Y() - xyz[1],2) + pow(d2.Z() - xyz[2],2) ); 

      float e1 = s1.E() ;
      float e2 = s2.E() ;

      _true_gamma_e_min = e1 < e2 ? e1 : e2 ;
      _true_gamma_e_max = e1 > e2 ? e1 : e2 ;

      _true_RL_minE = e1 < e2 ? radL1 : radL2 ;
      _true_RL_maxE = e1 > e2 ? radL1 : radL2 ;

      _true_asym = e1 > e2 ? e2/e1 : e1/e2 ;
      _true_angle = acos( dot / mag1 / mag2 ); 
    }

    if ( shr_ids.size() == 3 ) std::cout<<"energy: "<<_true_gamma_e_max <<std::endl ;
    if ( _true_gamma_e_max == -999 ) std::cout<<"shr_ids size: "<<shr_ids.size()<<std::endl ;

    // Also, get muon unfo
    auto ev_mctrk = storage->get_data<event_mctrack>("mcreco");
    if ( !ev_mctrk || ev_mctrk->size() == 0 ) {
      std::cout<<"No mctracks :( !"<<std::endl;
      return false;
    }

    for( int i = 0 ; i < ev_mctrk->size(); i++){
      auto mct = ev_mctrk->at(i);
      
      if ( mct.Origin() == 1 && mct.PdgCode() == 13 ){
        _true_mu_mom = mct.Start().Momentum().P() ;
        _true_mu_len = sqrt( pow(mct.End().X() - mct.Start().X(),2) +  
	                     pow(mct.End().Y() - mct.Start().Y(),2) +  
			     pow(mct.End().Z() - mct.Start().Z(),2) ) ;
        _true_mu_theta = mct.Start().Momentum().Theta() ;
        _true_mu_phi = mct.Start().Momentum().Phi() ;
      }
    }

   _pi0_tree->Fill();
      
    return true;
  }

  bool MCVariableStudy::finalize() {


  std::cout<<"Found N signals: "<<_n_signals<<std::endl;

    if(_fout) { 
      _fout->cd(); 
      _pi0_tree->Write(); 
      _gamma_tree->Write(); 
      }
  
    return true;
  }

}
#endif
