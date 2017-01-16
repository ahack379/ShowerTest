#ifndef LARLITE_CONTAINMENTSTUDY_CXX
#define LARLITE_CONTAINMENTSTUDY_CXX

#include "ContainmentStudy.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include <algorithm>

namespace larlite {

  bool ContainmentStudy::initialize() {
    
    if( !_pi0_tree ){
      _pi0_tree = new TTree("pi0_tree","");
      _pi0_tree->Branch("true_pi0_e",&_true_pi0_e,"true_pi0_e/F");
      _pi0_tree->Branch("true_angle",&_true_angle,"true_angle/F");
      _pi0_tree->Branch("true_asym",&_true_asym,"true_asym/F");
      _pi0_tree->Branch("reco_pi0_e",&_reco_pi0_e,"reco_pi0_e/F");
      _pi0_tree->Branch("event",&_event,"event/I");
      }

    if( !_gamma_tree ){
      _gamma_tree = new TTree("gamma_tree","");
      _gamma_tree->Branch("true_gamma_e",&_true_gamma_e,"true_gamma_e/F");
      _gamma_tree->Branch("reco_gamma_e",&_reco_gamma_e,"reco_gamma_e/F");
      _gamma_tree->Branch("true_reco_dot",&_true_reco_dot,"true_reco_dot/F");
      _gamma_tree->Branch("true_rad_l",&_true_rad_l,"true_rad_l/F");
      _gamma_tree->Branch("reco_rad_l",&_reco_rad_l,"reco_rad_l/F");
      }

    if( !_energy_tree ){
      _energy_tree = new TTree("energy_tree","");
      _energy_tree->Branch("containment",&_containment,"containment_e/F");
      _energy_tree->Branch("true_pi0_e",&_true_pi0_e,"true_pi0_e/F");
      _energy_tree->Branch("sum_gamma_e",&_sum_gamma_e,"sum_gamma_e/F");
      _energy_tree->Branch("mass",&_mass,"mass/F");
      _energy_tree->Branch("true_angle",&_true_angle,"true_angle/F");
      }

    _event = -1;

    return true;
  }

  void ContainmentStudy::Clear(){
    _true_pi0_e = -999;
    _true_angle = -9;
    _true_asym = -9;
    _reco_pi0_e = -999;

    _n_true_pi0 = 0;

    _true_gamma_e = -9;
    _reco_gamma_e = -9;
    _true_rad_l = -9;
    _reco_rad_l = -9;
    _true_reco_dot = -9;

    _containment = -1 ;
    _sum_gamma_e = 0;
    _mass = 0;
    }
  
  bool ContainmentStudy::analyze(storage_manager* storage) {

    _event ++; 
    Clear() ;

    std::cout<<"\nEvent is : "<<_event <<std::endl ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ){
      std::cout<<"No mctruth..." <<std::endl;
      return false;
      }

    auto parts = ev_mctruth->at(0).GetParticles();

    bool pi0 = false;
    bool mu  = false ;
    
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
    
    if ( _n_true_pi0 != 1 ) return false; 
      
   // Replace pi0 energy with combined shower energy 
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco"); 
    if(!ev_mcs || !ev_mcs->size() ){
      std::cout<<"No mcshower..." <<std::endl;
      return false;
      }

    std::vector<int> shr_ids;
     
    for ( int si = 0; si < ev_mcs->size(); si++){ 

      auto s = ev_mcs->at(si);

      if( s.MotherPdgCode() != 111 ) continue; 
      
      auto st = s.Start();
      auto dist = sqrt( pow(st.X() - start[0],2) + pow(st.Y() - start[1],2) + pow(st.Z() - start[2],2) );
      
      if ( dist < 0.001 ){
        _sum_gamma_e += s.DetProfile().E();
	shr_ids.emplace_back(si);
	}
      }

    _containment = _sum_gamma_e / _true_pi0_e ; 
    
    if( shr_ids.size() == 2 ){

      auto s1 = ev_mcs->at(shr_ids[0]).DetProfile();
      auto s2 = ev_mcs->at(shr_ids[1]).DetProfile();

      auto mag1 = sqrt( s1.Px()*s1.Px()+s1.Py()*s1.Py()+s1.Pz()*s1.Pz() );
      auto mag2 = sqrt( s2.Px()*s2.Px()+s2.Py()*s2.Py()+s2.Pz()*s2.Pz() );
      auto dot = s1.Px()*s2.Px() + s1.Py()*s2.Py() + s1.Pz()*s2.Pz() ;

      float e1 = s1.E() ;
      float e2 = s2.E() ;

      std::cout<<"DOT: "<<dot<<", "<<mag1<<", "<<mag2<<std::endl ; 

      if ( mag1 == 0 || mag2 == 0) return true; 
      //_true_asym = e1 > e2 ? e2/e1 : e1/e2 ;
      // true angle = cos (angle) 
      _true_angle = ( dot / mag1 / mag2 ); 
      _mass = sqrt(2 * e1 * e2 *(1. - _true_angle) );
      std::cout<<"Energy, angle, mass: "<<e1<<", "<<e2<<", "<<_true_angle<<", "<<_mass<<std::endl;
    }

    _energy_tree->Fill();
    //  
    //auto ev_s = storage->get_data<event_shower>("showerreco"); 
    //if( !ev_s || !ev_s->size() ){
    //  std::cout<<"No shower..." <<std::endl;
    //  return false;
    //  }

    //if( ev_s->size() < 2 ) return false ;
    //
    //float min_IP = 1e9;
    //int min_it = -1 ;

    //std::vector<std::pair<int,int>> CandidatePairs;
    //std::vector<int> cand_ids;

    ////std::cout<<ev_s->size()<<" Showers : "<<std::endl; 

    //for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

    //        auto const& shr1 = ev_s->at(s1);

    //        for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

    //            if (s2 <= s1) continue;

    //            auto const& shr2 = ev_s->at(s2);

    //            geoalgo::Vector_t dflipA(-1.*shr1.Direction()) ; 
    //            geoalgo::Vector_t dflipB(-1.*shr2.Direction()) ;

    //            // Make the backwards projection for the showers
    //            auto bksa = ::geoalgo::HalfLine_t(shr1.ShowerStart(),dflipA);
    //            auto bksb = ::geoalgo::HalfLine_t(shr2.ShowerStart(),dflipB);

    //            // Calc the Opening angle of the showers
    //            double Oangle = acos( shr1.Direction().Dot(shr2.Direction())) ; 

    //            // Calc the vertex point of the two showers. the true designated backwards project
    //            geoalgo::Point_t vertex(3);
    //            
    //            auto st1 = shr1.ShowerStart();
    //            auto st2 = shr2.ShowerStart();
    //            auto dir1 = shr1.Direction();
    //            auto dir2 = shr2.Direction();
    //            geoalgo::HalfLine_t shr1_hl(st1.X(),st1.Y(),st1.Z(),dir1.X(),dir1.Y(), dir1.Z() );
    //            geoalgo::HalfLine_t shr2_hl(st2.X(),st2.Y(),st2.Z(),dir2.X(),dir2.Y(), dir2.Z() );

    //            _geoAlgo.commonOrigin(shr1_hl, shr2_hl, vertex, true);

    //            // Calc Diretion of two correlated shower
    //            geoalgo::Vector_t momentum(3);// need to fill out
    //            geoalgo::Vector_t mom_vect(shr2.Direction()*shr1.Energy(2) +shr1.Direction()*shr2.Energy(2)) ;
    //            mom_vect.Normalize();
    //            momentum = mom_vect * sqrt(pow(shr1.Energy(),2)+pow(shr2.Energy(),2)+2*shr2.Energy()*shr1.Energy()*cos(Oangle));

    //            //===========================================
    //            auto _RIP = pow(_geoAlgo.SqDist(bksa,bksb),0.5);
    //            auto _RadL_A = vertex.Dist(shr1.ShowerStart());
    //            auto _RadL_B = vertex.Dist(shr2.ShowerStart());

    //            if(Oangle < 0.35){
    //              //std::cout<<"Bad Angle: "<< Oangle<<std::endl ; 
    //              continue;
    //             }

    //            if(pow(_geoAlgo.SqDist(bksa,bksb),0.5) > 4){
    //               //std::cout<<"Ip? "<<pow(_geoAlgo.SqDist(bksa,bksb),0.5)<<std::endl ; 
    //               continue;
    //              }

    //            if(_RadL_A > 62. || _RadL_B > 62. ){ 
    //               //std::cout<<"Rad Length : "<<_RadL_A<<", "<<_RadL_B<<std::endl ;
    //               continue; 
    //                }

    //            if( vertex[0] > 256.35 || vertex[0] < 0. || vertex[1] > 116.5  
    //              || vertex[1] < -116.5 || vertex[2] > 1036.8 || vertex[2] < 0. ){

    //              std::cout<<"Vertex out of bounds : "<<vertex[0]<<", "<<vertex[1]<<", "<<vertex[2]<<std::endl ;
    //              continue;
    //              }

    //            // Bunch of cuts
    //            if( _RIP < min_IP ){
    //              min_IP = _RIP ;
    //              min_it = CandidatePairs.size();
    //              }

    //            CandidatePairs.push_back(std::make_pair(s1,s2)); 
    //    	cand_ids.emplace_back(s1);
    //    	cand_ids.emplace_back(s2);

    //        }// shower ID 2 
    //    }// shower ID 1 

    //   std::cout<<"Cand size: "<<cand_ids.size()<<std::endl ;
    //   // Fill trees 
    //   if( cand_ids.size() == 2 && min_it != -1){ 
    //
    //       auto s1 = ev_s->at(cand_ids[0]) ;
    //       auto s2 = ev_s->at(cand_ids[1]) ;
    //      
    //       _reco_pi0_e = s1.Energy(2) + s2.Energy(2); 
    //       _pi0_tree->Fill();

    //       std::vector<std::pair<int,int>> mc_to_reco ;
    //       float max_dot = -1e9;
    //       int max_mcs(-1), max_recos(-1), min_mcs(-1), min_recos(-1);

    //       // Match showers
    //       for( auto const & mc_id : shr_ids ){
    //         auto mcs_i = ev_mcs->at(mc_id);
    //         auto mag_mcs = sqrt( pow(mcs_i.Start().Px(),2) + pow(mcs_i.Start().Py(),2) + pow(mcs_i.Start().Pz(),2) ); 
    //         
    //         for( auto const & reco_id : cand_ids ){

    //           auto recos_i = ev_s->at(reco_id) ;

    //           auto mag_reco = sqrt( pow(recos_i.Direction().Px(),2) + pow(recos_i.Direction().Py(),2) + pow(recos_i.Direction().Pz(),2) ); 
    //           auto dot = mcs_i.Start().Px() * recos_i.Direction().Px() + 
    //                      mcs_i.Start().Py() * recos_i.Direction().Py() + 
    //    		  mcs_i.Start().Pz() * recos_i.Direction().Pz() ;
    //           dot /= ( mag_mcs * mag_reco );

    //           if ( dot > max_dot){
    //             max_dot = dot;
    //             max_mcs = mc_id;
    //             max_recos = reco_id;
    //            }
    //          } 
    //        }
    //        
    //       min_mcs = ( max_mcs == shr_ids[0] ? shr_ids[1] : shr_ids[0] ) ;
    //       min_recos = ( max_recos == cand_ids[0] ? cand_ids[1] : cand_ids[0] ) ;

    //       auto mcs_1 = ev_mcs->at(max_mcs) ;
    //       auto recos_1 = ev_s->at(max_recos) ;

    //       _true_gamma_e = mcs_1.Start().E();
    //       _reco_gamma_e = recos_1.Energy(2);
    //       _true_reco_dot = max_dot;
    //       _gamma_tree->Fill();

    //       auto mcs_2 = ev_mcs->at(min_mcs) ;
    //       auto recos_2 = ev_s->at(min_recos) ;
    //       auto dot = mcs_2.Start().Px() * recos_2.Direction().Px() + mcs_2.Start().Py() * recos_2.Direction().Py() + mcs_2.Start().Pz() * recos_2.Direction().Pz() ;
    //       auto mag_reco = sqrt( pow(recos_2.Direction().Px(),2) + pow(recos_2.Direction().Py(),2) + pow(recos_2.Direction().Pz(),2) ); 
    //       auto mag_mcs = sqrt( pow(mcs_2.Start().Px(),2) + pow(mcs_2.Start().Py(),2) + pow(mcs_2.Start().Pz(),2) ); 
    //       dot /= (mag_mcs * mag_reco);

    //       _true_gamma_e = mcs_2.Start().E();
    //       _reco_gamma_e = recos_2.Energy(2);
    //       _true_reco_dot = dot;
    //       _gamma_tree->Fill();
    //     }
  
    return true;
  }

  bool ContainmentStudy::finalize() {

    if(_fout) { 
      _fout->cd(); 
      //_pi0_tree->Write(); 
      //_gamma_tree->Write(); 
      _energy_tree->Write();
      }
  
    return true;
  }

}
#endif
