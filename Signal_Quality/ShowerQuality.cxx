#ifndef LARLITE_SHOWERQUALITY_CXX
#define LARLITE_SHOWERQUALITY_CXX

#include "ShowerQuality.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/vertex.h"
#include "DataFormat/track.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/shower.h"
#include "DataFormat/mcshower.h"
#include <typeinfo> 

// This module checks the final selected sample of CCpi0's for true signal CCpi0 events
// For each of these events, we will study the quality of the pieces:
// vertex resolution, track start point + directional resolution, start 
// point and energy resolution of each reconstruction shower

namespace larlite {

  bool ShowerQuality::initialize() {

    _pi0s = 0;
    _fv = 0;

    if(!_vtx_tree){
     _vtx_tree = new TTree("vtx_tree","vtx_tree"); 
     _vtx_tree->Branch("event",&_event,"event/I"); 
     _vtx_tree->Branch("vtx_diff",&_vtx_diff,"vtx_diff/F"); 
     _vtx_tree->Branch("mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/F"); 
     _vtx_tree->Branch("mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/F"); 
     _vtx_tree->Branch("mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/F"); 
     _vtx_tree->Branch("reco_vtx_x",&_reco_vtx_x,"reco_vtx_x/F"); 
     _vtx_tree->Branch("reco_vtx_y",&_reco_vtx_y,"reco_vtx_y/F"); 
     _vtx_tree->Branch("reco_vtx_z",&_reco_vtx_z,"reco_vtx_z/F"); 

     _vtx_tree->Branch("reco_shr1_e",&_reco_shr1_e,"reco_shr1_e/F"); 
     _vtx_tree->Branch("mc_shr1_e",&_mc_shr1_e,"mc_shr1_e/F"); 
     _vtx_tree->Branch("shr1_dot",&_shr1_dot,"shr1_dot/F"); 
     _vtx_tree->Branch("reco_shr1_st_x",&_reco_shr1_st_x,"reco_shr1_st_x/F"); 
     _vtx_tree->Branch("reco_shr1_st_y",&_reco_shr1_st_y,"reco_shr1_st_y/F"); 
     _vtx_tree->Branch("reco_shr1_st_z",&_reco_shr1_st_z,"reco_shr1_st_z/F"); 
     _vtx_tree->Branch("mc_shr1_st_x",&_mc_shr1_st_x,"mc_shr1_st_x/D"); 
     _vtx_tree->Branch("mc_shr1_st_y",&_mc_shr1_st_y,"mc_shr1_st_y/D"); 
     _vtx_tree->Branch("mc_shr1_st_z",&_mc_shr1_st_z,"mc_shr1_st_z/D"); 
     _vtx_tree->Branch("shr1_st_diff",&_shr1_st_diff,"shr1_st_diff/F"); 

     _vtx_tree->Branch("reco_shr2_e",&_reco_shr2_e,"reco_shr2_e/F"); 
     _vtx_tree->Branch("mc_shr2_e",&_mc_shr2_e,"mc_shr2_e/F"); 
     _vtx_tree->Branch("shr2_dot",&_shr2_dot,"shr2_dot/F"); 
     _vtx_tree->Branch("reco_shr2_st_x",&_reco_shr2_st_x,"reco_shr2_st_x/F"); 
     _vtx_tree->Branch("reco_shr2_st_y",&_reco_shr2_st_y,"reco_shr2_st_y/F"); 
     _vtx_tree->Branch("reco_shr2_st_z",&_reco_shr2_st_z,"reco_shr2_st_z/F"); 
     _vtx_tree->Branch("mc_shr2_st_x",&_mc_shr2_st_x,"mc_shr2_st_x/D"); 
     _vtx_tree->Branch("mc_shr2_st_y",&_mc_shr2_st_y,"mc_shr2_st_y/D"); 
     _vtx_tree->Branch("mc_shr2_st_z",&_mc_shr2_st_z,"mc_shr2_st_z/D"); 
     _vtx_tree->Branch("shr2_st_diff",&_shr2_st_diff,"shr2_st_diff/F"); 

    }   

    SetXOffset(0.0);
    //_offset = 0.;
    _event = 0;
    _bad_events = 0;

    _SCE = new larutil::SpaceChargeMicroBooNE();
    _time2cm = larutil::GeometryHelper::GetME()->TimeToCm();

    return true;
  }
  
  bool ShowerQuality::analyze(storage_manager* storage) {

   // Iterate event number here to make event display checking convenient
   _event++;

   // Make sure no duplicate events are considered; if encounter the 2nd
   // event in a duplicate pair, do not consider the event for tree filling
   auto it = _map.find(storage->subrun_id());
   bool foundit = false;
   if( it != _map.end() ){
    while ( it->first == storage->subrun_id() ){  
      auto temp_event = it->second ; 
      if( temp_event == storage->event_id() )
        foundit = true;

      it++; 
      }   
     if ( !foundit)
      _map.emplace(storage->subrun_id(), storage->event_id() );

     else return false ;
    }   
   else 
      _map.emplace(storage->subrun_id(), storage->event_id() );


   std::cout<<"\n\nON EVENT : "<<_event - 1 <<std::endl ;

   auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) { 
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
      }   

   auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex"); 
    if(!ev_vtx || !ev_vtx->size() ) { 
      std::cout<<"Event has no vertex info "<<std::endl;
      return false;
      }

   auto ev_mcshr = storage->get_data<event_mcshower>("mcreco"); 
    if(!ev_mcshr || !ev_mcshr->size() ) { 
      std::cout<<"Event has no mcshower info "<<std::endl;
      return false;
      }

   auto ev_s = storage->get_data<event_shower>("showerreco");
   if( !ev_s || !ev_s->size() ){
     std::cout<<"Event has no shower info..." <<std::endl;
     return false;
     }


    // First get truth information so we can select only true signal from this final sample
    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    auto traj = nu.Nu().Trajectory();
    auto xvtx = traj.at(traj.size() - 1).X();
    auto yvtx = traj.at(traj.size() - 1).Y();
    auto zvtx = traj.at(traj.size() - 1).Z();
    auto e = traj.at(traj.size() - 1).E();

    // Only want to consider true ccpi0 events with vertices in the FV
    if( xvtx < 20 || xvtx > 236.35 || yvtx > 96.5 || yvtx < -96.5 || zvtx < 10 || zvtx > 1026.8 ){
        std::cout<<"EVENT : "<<_event -1 <<" is OUT OF FV "<<std::endl ;
	_fv++;
        return false;
	}

    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0;
    int n_mu = 0;
    
    for ( auto const & p : parts ){
    
        if( p.StatusCode() == 1 && p.PdgCode() == 111 )
          n_pi0 ++; 
        if( p.StatusCode() == 1 && p.PdgCode() == 13 )
          n_mu ++; 
     }   
    
    
    if ( n_pi0 == 1 && n_mu == 1 && e > 0.5){ 

      // If we get here, we found a signal event! Now start storing shit 

      // Correct for space charge effect
      auto tvtx = traj.at(traj.size() - 1).T(); // ns
      auto vtxtick = (tvtx / 1000.) * 2.; // time in tick :
      auto vtxtimecm = vtxtick * _time2cm; // time in cm :
      // Correct for space charge effects in simulation
      auto sce_corr = _SCE->GetPosOffsets(xvtx,yvtx,zvtx);

      // Store first MC-Reco vertex comparisons
      _mc_vtx_x = xvtx + vtxtimecm + _offset - sce_corr.at(0); 
      _mc_vtx_y = yvtx + sce_corr.at(1); 
      _mc_vtx_z = zvtx + sce_corr.at(2);

      auto reco_vtx = ev_vtx->at(0);
      _reco_vtx_x = reco_vtx.X(); 
      _reco_vtx_y = reco_vtx.Y(); 
      _reco_vtx_z = reco_vtx.Z(); 

      _vtx_diff = sqrt(pow(_reco_vtx_x - _mc_vtx_x, 2) + pow(_reco_vtx_y - _mc_vtx_y, 2) 
                     + pow(_reco_vtx_z - _mc_vtx_z, 2) ); 


      // Now deal with shower comparisons
      std::vector<int> shr_ids;

      for ( int si = 0; si < ev_mcshr->size(); si++){

        auto s = ev_mcshr->at(si);
        auto st = s.Start();
        auto dist = sqrt( pow(st.X() - xvtx,2) + pow(st.Y() - yvtx,2) + pow(st.Z() - zvtx,2) );

        if ( dist < 0.01 && s.DetProfile().E() > 0 && s.MotherPdgCode() == 111 ){
	  std::cout<<"ENERGY : "<<s.DetProfile().E()<<std::endl ;
	  std::cout<<"START INFO : "<<s.DetProfile().X() <<", "<<s.DetProfile().Y()<<", "<<s.DetProfile().Z()<<std::endl ;
	  _pi0s ++; 
          shr_ids.emplace_back(si) ;
          }
      }

     if (shr_ids.size() != 2 ) { //< 2){
       std::cout<<"N SHOWERS! "<<shr_ids.size()<<", Event: "<<_event-1<<std::endl ;
       return false;
     }

    std::vector<std::pair<int,int>> CandidatePairs;
    std::vector<int> cand_ids;

    for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

            auto const& shr1 = ev_s->at(s1);

            for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

                if (s2 <= s1) continue;

                auto const& shr2 = ev_s->at(s2);

                geoalgo::Vector_t dflipA(-1.*shr1.Direction()) ;
                geoalgo::Vector_t dflipB(-1.*shr2.Direction()) ;

                // Make the backwards projection for the showers
                auto bksa = ::geoalgo::HalfLine_t(shr1.ShowerStart(),dflipA);
                auto bksb = ::geoalgo::HalfLine_t(shr2.ShowerStart(),dflipB);

                // Calc the Opening angle of the showers
                double Oangle = acos( shr1.Direction().Dot(shr2.Direction())) ;

                // Calc the vertex point of the two showers. the true designated backwards project
                geoalgo::Point_t vertex(3);

                auto st1 = shr1.ShowerStart();
                auto st2 = shr2.ShowerStart();
                auto dir1 = shr1.Direction();
                auto dir2 = shr2.Direction();
                geoalgo::HalfLine_t shr1_hl(st1.X(),st1.Y(),st1.Z(),dir1.X(),dir1.Y(), dir1.Z() );
                geoalgo::HalfLine_t shr2_hl(st2.X(),st2.Y(),st2.Z(),dir2.X(),dir2.Y(), dir2.Z() );

                _geoAlgo.commonOrigin(shr1_hl, shr2_hl, vertex, true);

                // Calc Diretion of two correlated shower
                geoalgo::Vector_t momentum(3);// need to fill out
                geoalgo::Vector_t mom_vect(shr2.Direction()*shr1.Energy(2) +shr1.Direction()*shr2.Energy(2)) ;
                mom_vect.Normalize();
                momentum = mom_vect * sqrt(pow(shr1.Energy(),2)+pow(shr2.Energy(),2)+2*shr2.Energy()*shr1.Energy()*cos(Oangle));

                //===========================================
                auto _RadL_A = vertex.Dist(shr1.ShowerStart());
                auto _RadL_B = vertex.Dist(shr2.ShowerStart());

                if(Oangle < 0.35){ continue; }
                if(pow(_geoAlgo.SqDist(bksa,bksb),0.5) > 4) continue;
                if(_RadL_A > 62. || _RadL_B > 62. ){ continue; }

                if( vertex[0] > 256.35 || vertex[0] < 0. || vertex[1] > 116.5
                  || vertex[1] < -116.5 || vertex[2] > 1036.8 || vertex[2] < 0. ){ continue; }

                CandidatePairs.push_back(std::make_pair(s1,s2));
                cand_ids.emplace_back(s1);
                cand_ids.emplace_back(s2);

            }// shower ID 2 
        }// shower ID 1 


       if( cand_ids.size() == 2 && shr_ids.size() > 1){

           auto s1 = ev_s->at(cand_ids[0]) ;
           auto s2 = ev_s->at(cand_ids[1]) ;

           std::vector<std::pair<int,int>> mc_to_reco ;
           float max_dot = -1e9;
           int max_mcs(-1), max_recos(-1), min_mcs(-1), min_recos(-1);

           // Match showers
           for( auto const & mc_id : shr_ids ){
             auto mcs_i = ev_mcshr->at(mc_id);
             auto mag_mcs = sqrt( pow(mcs_i.DetProfile().Px(),2) + pow(mcs_i.DetProfile().Py(),2) + pow(mcs_i.DetProfile().Pz(),2) );

             for( auto const & reco_id : cand_ids ){

               auto recos_i = ev_s->at(reco_id) ;

               auto mag_reco = sqrt( pow(recos_i.Direction().Px(),2) + pow(recos_i.Direction().Py(),2) + pow(recos_i.Direction().Pz(),2) );
               auto dot = mcs_i.DetProfile().Px() * recos_i.Direction().Px() +
                          mcs_i.DetProfile().Py() * recos_i.Direction().Py() +
                          mcs_i.DetProfile().Pz() * recos_i.Direction().Pz() ;
               dot /= ( mag_mcs * mag_reco );
	       if ( fabs(dot) > 1 ) 
	       std::cout<<"DOT ! " <<dot <<std::endl ;

               if ( dot > max_dot){
                 max_dot = dot;
                 max_mcs = mc_id;
                 max_recos = reco_id;
                }
              }
            }

           min_mcs = ( max_mcs == shr_ids[0] ? shr_ids[1] : shr_ids[0] ) ;
           min_recos = ( max_recos == cand_ids[0] ? cand_ids[1] : cand_ids[0] ) ;

           auto mcs_1 = ev_mcshr->at(max_mcs) ;
           auto recos_1 = ev_s->at(max_recos) ;

           _reco_shr1_e= recos_1.Energy(2);
           _mc_shr1_e  = mcs_1.DetProfile().E();
           _shr1_dot   = max_dot;

           _reco_shr1_st_x = recos_1.ShowerStart().X();
           _reco_shr1_st_y = recos_1.ShowerStart().Y();
           _reco_shr1_st_z = recos_1.ShowerStart().Z();

           auto st_x = mcs_1.DetProfile().X();
           auto st_y = mcs_1.DetProfile().Y();
           auto st_z = mcs_1.DetProfile().Z();

           auto sce_corr = _SCE->GetPosOffsets(st_x,st_y,st_z);

           _mc_shr1_st_x = st_x + vtxtimecm + _offset - sce_corr.at(0); 
           _mc_shr1_st_y = st_y + sce_corr.at(1); 
           _mc_shr1_st_z = st_z + sce_corr.at(2);
 
           _shr1_st_diff = sqrt( pow(_reco_shr1_st_x - _mc_shr1_st_x,2) + pow(_reco_shr1_st_y - _mc_shr1_st_y,2) + 
                                  pow(_reco_shr1_st_z - _mc_shr1_st_z,2) ); 

           auto mcs_2 = ev_mcshr->at(min_mcs) ;
           auto recos_2 = ev_s->at(min_recos) ;
           auto dot = mcs_2.Start().Px() * recos_2.Direction().Px() + mcs_2.Start().Py() * recos_2.Direction().Py() + mcs_2.Start().Pz() * recos_2.Direction().Pz() ;
           auto mag_reco = sqrt( pow(recos_2.Direction().Px(),2) + pow(recos_2.Direction().Py(),2) + pow(recos_2.Direction().Pz(),2) );
           auto mag_mcs = sqrt( pow(mcs_2.Start().Px(),2) + pow(mcs_2.Start().Py(),2) + pow(mcs_2.Start().Pz(),2) );
           dot /= (mag_mcs * mag_reco);

           _reco_shr2_e= recos_2.Energy(2);
           _mc_shr2_e  = mcs_2.DetProfile().E();
           _shr2_dot   = dot;

           _reco_shr2_st_x = recos_2.ShowerStart().X();
           _reco_shr2_st_y = recos_2.ShowerStart().Y();
           _reco_shr2_st_z = recos_2.ShowerStart().Z();

           auto st2_x = mcs_2.DetProfile().X();
           auto st2_y = mcs_2.DetProfile().Y();
           auto st2_z = mcs_2.DetProfile().Z();

           sce_corr = _SCE->GetPosOffsets(st2_x,st2_y,st2_z);

           _mc_shr2_st_x = st2_x + vtxtimecm + _offset - sce_corr.at(0); 
           _mc_shr2_st_y = st2_y + sce_corr.at(1); 
           _mc_shr2_st_z = st2_z + sce_corr.at(2);
           
           _shr2_st_diff = sqrt( pow(_reco_shr2_st_x - _mc_shr2_st_x,2) + pow(_reco_shr2_st_y - _mc_shr2_st_y,2) + 
                                  pow(_reco_shr2_st_z - _mc_shr2_st_z,2) ); 

	   if ( _shr1_st_diff > 8 || _shr2_st_diff > 8 ){ 
	     std::cout<<"FAILED EVENT: "<<_event-1<<", "<<_shr1_st_diff<<", "<<_shr2_st_diff<<std::endl ;
	     _bad_events++;
	     }
           else
             std::cout<<"SUCCESS! "<<_event -1 <<", "<<shr_ids.size()<<std::endl ;

           _vtx_tree->Fill();

         }

      _event_list.emplace_back(_event-1);

    }
    else{
      std::cout<<"Not a CCPI0 EVENT ! " <<std::endl ;
      }

    return true;
  }

  bool ShowerQuality::finalize() {

   std::cout<<"Event list: "<<_event_list.size()<<std::endl ; 
   std::cout<<"Bad events: "<<_bad_events<<std::endl ;
   std::cout<<"Out of FV: "<<_fv<<std::endl ;

   std::cout<<"N showers with pi0s : (should be 2x evnetlist) "<<_pi0s<<std::endl ;

   if(_fout){
     _fout->cd();
     _vtx_tree->Write();
   }
  
    return true;
  }

}
#endif
