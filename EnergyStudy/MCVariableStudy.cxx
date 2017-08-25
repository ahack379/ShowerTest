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

    return true;
  }

  void MCVariableStudy::Clear(){
    _true_pi0_e = -999;
    _true_angle = -9;
    _true_asym = -9;
    _reco_pi0_e = -999;
    _true_gamma_e_max = 0 ;
    _true_gamma_e_min = 0 ;
    _true_pi0_mom = -999;
    _true_nu_e = -999;
    _true_RL_maxE = 0;
    _true_RL_minE = 0;

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

    std::cout<<"\nEvent is : "<<_event <<std::endl ;

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

    if( xyz[0] < 20 || xyz[0] > 236.25 || xyz[1] < -96.5 || xyz[1] > 96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
      return false ;

    auto parts = ev_mctruth->at(0).GetParticles();

    int n_mes = 0;
    int n_lep = 0;
    int n_mu = 0;
    
    std::vector<float> start(3,0) ;
    
    for ( auto const & p : parts ){
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){
        _n_true_pi0 ++;
        _true_pi0_e = p.Trajectory().at(0).E()*1000; 
        start[0] = p.Trajectory().at(0).X();
        start[1] = p.Trajectory().at(0).Y();
        start[2] = p.Trajectory().at(0).Z();
	_true_pi0_mom = sqrt( pow(p.Trajectory().at(0).Px(),2) + pow(p.Trajectory().at(0).Py(),2) + pow(p.Trajectory().at(0).Pz(),2) )*1000; 
      }

      if( p.StatusCode() == 1 && p.PdgCode() == 13 )
        n_mu ++; 

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 211 || abs(p.PdgCode()) == 321 ||  p.PdgCode() == 130 
                                   || p.PdgCode() == 310 || abs(p.PdgCode()) == 311 ) ){
        n_mes ++; 
      }   

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 11 || p.PdgCode() == -13))
        n_lep ++; 

    }

    if( n_mu != 1 || _n_true_pi0 != 1) return false ; // || n_lep != 0 || n_mes != 0) return false;

    _true_nu_e = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().at(0).E() ;
      
   // Replace pi0 energy with combined shower energy 
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco"); 
    if(!ev_mcs || !ev_mcs->size() ){
      std::cout<<"No mcshower..." <<std::endl;
      return false;
    }

    std::vector<int> shr_ids;
     
    for ( int si = 0; si < ev_mcs->size(); si++){ 

      auto s = ev_mcs->at(si);

      if( s.PdgCode() != 22 ) continue; //|| abs(s.MotherPdgCode()) == 13 ) continue;
      
      auto st = s.Start();
      auto dist = sqrt( pow(st.X() - start[0],2) + pow(st.Y() - start[1],2) + pow(st.Z() - start[2],2) );
      
      if ( dist < 0.001 ){
        shr_ids.emplace_back(si) ;
        //std::cout<<"Dist: "<<dist<<", energies shr: "<<ev_mcs->at(si).Start().E()<<std::endl;
	}
    }
    
    if( shr_ids.size() == 2 ){

      auto s1 = ev_mcs->at(shr_ids[0]).Start();
      auto s2 = ev_mcs->at(shr_ids[1]).Start();

      auto d1 = ev_mcs->at(shr_ids[0]).DetProfile();
      auto d2 = ev_mcs->at(shr_ids[1]).DetProfile();

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
      
   // auto ev_s = storage->get_data<event_shower>("showerreco"); 
   // if( !ev_s || !ev_s->size() ){
   //   std::cout<<"No shower..." <<std::endl;
   //   return false;
   //   }

   // if( ev_s->size() < 2 ) return false ;
   // 
   // float min_IP = 1e9;
   // int min_it = -1 ;

   // std::vector<std::pair<int,int>> CandidatePairs;
   // std::vector<int> cand_ids;

    //std::cout<<ev_s->size()<<" Showers : "<<std::endl; 

   // for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

   //         auto const& shr1 = ev_s->at(s1);

   //         for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

   //             if (s2 <= s1) continue;

   //             auto const& shr2 = ev_s->at(s2);

   //             geoalgo::Vector_t dflipA(-1.*shr1.Direction()) ; 
   //             geoalgo::Vector_t dflipB(-1.*shr2.Direction()) ;

   //             // Make the backwards projection for the showers
   //             auto bksa = ::geoalgo::HalfLine_t(shr1.ShowerStart(),dflipA);
   //             auto bksb = ::geoalgo::HalfLine_t(shr2.ShowerStart(),dflipB);

   //             // Calc the Opening angle of the showers
   //             double Oangle = acos( shr1.Direction().Dot(shr2.Direction())) ; 

   //             // Calc the vertex point of the two showers. the true designated backwards project
   //             geoalgo::Point_t vertex(3);
   //             
   //             auto st1 = shr1.ShowerStart();
   //             auto st2 = shr2.ShowerStart();
   //             auto dir1 = shr1.Direction();
   //             auto dir2 = shr2.Direction();
   //             geoalgo::HalfLine_t shr1_hl(st1.X(),st1.Y(),st1.Z(),dir1.X(),dir1.Y(), dir1.Z() );
   //             geoalgo::HalfLine_t shr2_hl(st2.X(),st2.Y(),st2.Z(),dir2.X(),dir2.Y(), dir2.Z() );

   //             _geoAlgo.commonOrigin(shr1_hl, shr2_hl, vertex, true);

   //             // Calc Diretion of two correlated shower
   //             geoalgo::Vector_t momentum(3);// need to fill out
   //             geoalgo::Vector_t mom_vect(shr2.Direction()*shr1.Energy(2) +shr1.Direction()*shr2.Energy(2)) ;
   //             mom_vect.Normalize();
   //             momentum = mom_vect * sqrt(pow(shr1.Energy(),2)+pow(shr2.Energy(),2)+2*shr2.Energy()*shr1.Energy()*cos(Oangle));

   //             //===========================================
   //             auto _RIP = pow(_geoAlgo.SqDist(bksa,bksb),0.5);
   //             auto _RadL_A = vertex.Dist(shr1.ShowerStart());
   //             auto _RadL_B = vertex.Dist(shr2.ShowerStart());

   //             if(Oangle < 0.35){
   //               //std::cout<<"Bad Angle: "<< Oangle<<std::endl ; 
   //               continue;
   //              }

   //             if(pow(_geoAlgo.SqDist(bksa,bksb),0.5) > 4){
   //                //std::cout<<"Ip? "<<pow(_geoAlgo.SqDist(bksa,bksb),0.5)<<std::endl ; 
   //                continue;
   //               }

   //             if(_RadL_A > 62. || _RadL_B > 62. ){ 
   //                //std::cout<<"Rad Length : "<<_RadL_A<<", "<<_RadL_B<<std::endl ;
   //                continue; 
   //                 }

   //             if( vertex[0] > 256.35 || vertex[0] < 0. || vertex[1] > 116.5  
   //               || vertex[1] < -116.5 || vertex[2] > 1036.8 || vertex[2] < 0. ){

   //               std::cout<<"Vertex out of bounds : "<<vertex[0]<<", "<<vertex[1]<<", "<<vertex[2]<<std::endl ;
   //               continue;
   //               }

   //             // Bunch of cuts
   //             if( _RIP < min_IP ){
   //               min_IP = _RIP ;
   //               min_it = CandidatePairs.size();
   //               }

   //             CandidatePairs.push_back(std::make_pair(s1,s2)); 
   //     	cand_ids.emplace_back(s1);
   //     	cand_ids.emplace_back(s2);

   //         }// shower ID 2 
   //     }// shower ID 1 

   //    std::cout<<"Cand size: "<<cand_ids.size()<<std::endl ;
   //    // Fill trees 
   //    if( cand_ids.size() == 2 && min_it != -1){ 
   // 
   //        auto s1 = ev_s->at(cand_ids[0]) ;
   //        auto s2 = ev_s->at(cand_ids[1]) ;
   //       
   //        _reco_pi0_e = s1.Energy(2) + s2.Energy(2); 
   //        _pi0_tree->Fill();

   //        std::vector<std::pair<int,int>> mc_to_reco ;
   //        float max_dot = -1e9;
   //        int max_mcs(-1), max_recos(-1), min_mcs(-1), min_recos(-1);

   //        // Match showers
   //        for( auto const & mc_id : shr_ids ){
   //          auto mcs_i = ev_mcs->at(mc_id);
   //          auto mag_mcs = sqrt( pow(mcs_i.Start().Px(),2) + pow(mcs_i.Start().Py(),2) + pow(mcs_i.Start().Pz(),2) ); 
   //          
   //          for( auto const & reco_id : cand_ids ){

   //            auto recos_i = ev_s->at(reco_id) ;

   //            auto mag_reco = sqrt( pow(recos_i.Direction().Px(),2) + pow(recos_i.Direction().Py(),2) + pow(recos_i.Direction().Pz(),2) ); 
   //            auto dot = mcs_i.Start().Px() * recos_i.Direction().Px() + 
   //                       mcs_i.Start().Py() * recos_i.Direction().Py() + 
   //     		  mcs_i.Start().Pz() * recos_i.Direction().Pz() ;
   //            dot /= ( mag_mcs * mag_reco );

   //            if ( dot > max_dot){
   //              max_dot = dot;
   //              max_mcs = mc_id;
   //              max_recos = reco_id;
   //             }
   //           } 
   //         }
   //         
   //        min_mcs = ( max_mcs == shr_ids[0] ? shr_ids[1] : shr_ids[0] ) ;
   //        min_recos = ( max_recos == cand_ids[0] ? cand_ids[1] : cand_ids[0] ) ;

   //        auto mcs_1 = ev_mcs->at(max_mcs) ;
   //        auto recos_1 = ev_s->at(max_recos) ;

   //        _true_gamma_e = mcs_1.Start().E();
   //        _reco_gamma_e = recos_1.Energy(2);
   //        _true_reco_dot = max_dot;
   //        _gamma_tree->Fill();
   //        std::cout<<"Reco and tru gamma e: "<<_reco_gamma_e<<", "<<_true_gamma_e<<std::endl ;
   //        std::cout<<"Start of shower: "<<mcs_1.Start().X()<<", "<<mcs_1.Start().Y() <<", "<<mcs_1.Start().Z()<<std::endl ;
   //        std::cout<<"End of shower: "<<mcs_1.End().X()<<", "<<mcs_1.End().Y() <<", "<<mcs_1.End().Z()<<std::endl ;

   //        auto mcs_2 = ev_mcs->at(min_mcs) ;
   //        auto recos_2 = ev_s->at(min_recos) ;
   //        auto dot = mcs_2.Start().Px() * recos_2.Direction().Px() + mcs_2.Start().Py() * recos_2.Direction().Py() + mcs_2.Start().Pz() * recos_2.Direction().Pz() ;
   //        auto mag_reco = sqrt( pow(recos_2.Direction().Px(),2) + pow(recos_2.Direction().Py(),2) + pow(recos_2.Direction().Pz(),2) ); 
   //        auto mag_mcs = sqrt( pow(mcs_2.Start().Px(),2) + pow(mcs_2.Start().Py(),2) + pow(mcs_2.Start().Pz(),2) ); 
   //        dot /= (mag_mcs * mag_reco);

   //        _true_gamma_e = mcs_2.Start().E();
   //        _reco_gamma_e = recos_2.Energy(2);
   //        _true_reco_dot = dot;
   //        _gamma_tree->Fill();

   //        std::cout<<"\nReco and tru gamma e: "<<_reco_gamma_e<<", "<<_true_gamma_e<<std::endl ;
   //        std::cout<<"Start of shower: "<<mcs_2.Start().X()<<", "<<mcs_2.Start().Y() <<", "<<mcs_2.Start().Z()<<std::endl ;
   //        std::cout<<"End of shower: "<<mcs_2.End().X()<<", "<<mcs_2.End().Y() <<", "<<mcs_2.End().Z()<<std::endl ;
   //      }
  
    return true;
  }

  bool MCVariableStudy::finalize() {

    if(_fout) { 
      _fout->cd(); 
      _pi0_tree->Write(); 
      _gamma_tree->Write(); 
      }
  
    return true;
  }

}
#endif
