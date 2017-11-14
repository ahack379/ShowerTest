#ifndef LARLITE_ACCEPTANCESTUDY_CXX
#define LARLITE_ACCEPTANCESTUDY_CXX

#include "AcceptanceStudy.h"
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

  bool AcceptanceStudy::initialize() {

    _infv_ccpi0 = 0;
    _thresh = 0;
    _dalitz = 0;
    _event = 0;
    _out_of_vol = 0;

    return true;
  }
  
  bool AcceptanceStudy::analyze(storage_manager* storage) {

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


   auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) { 
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
      }   

   auto ev_mcshr = storage->get_data<event_mcshower>("mcreco"); 
    if(!ev_mcshr || !ev_mcshr->size() ) { 
      std::cout<<"Event has no mcshower info "<<std::endl;
      return false;
      }

   //std::cout<<"\n\nON EVENT : "<<_event - 1 <<std::endl ;

    // First get truth information so we can select only true signal from this final sample
    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    auto traj = nu.Nu().Trajectory();
    auto xvtx = traj.at(traj.size() - 1).X();
    auto yvtx = traj.at(traj.size() - 1).Y();
    auto zvtx = traj.at(traj.size() - 1).Z();
    auto e = traj.at(traj.size() - 1).E();

    bool infv = true;

    // Only want to consider true ccpi0 events with vertices in the FV
    if( xvtx < 20 || xvtx > 236.35 || yvtx > 96.5 || yvtx < -96.5 || zvtx < 10 || zvtx > 1026.8 )
      infv = false;

    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0;
    int n_mu = 0;
    
    for ( auto const & p : parts ){
    
        if( p.StatusCode() == 1 && p.PdgCode() == 111 )
          n_pi0 ++; 
        if( p.StatusCode() == 1 && p.PdgCode() == 13 )
          n_mu ++; 
     }   
    
    if ( n_pi0 == 1 && n_mu == 1 && infv ){ 

      // If we get here, we found a signal event! 
      _infv_ccpi0++ ;

      std::vector<int> shr_ids;
      bool out = false;

      for ( int si = 0; si < ev_mcshr->size(); si++){

        auto s = ev_mcshr->at(si);
        auto st = s.Start();
        auto end = s.End();
        auto det = s.DetProfile();
        auto dist = sqrt( pow(st.X() - xvtx,2) + pow(st.Y() - yvtx,2) + pow(st.Z() - zvtx,2) );

        if ( dist < 0.0001 && s.MotherPdgCode() == 111 ){

         if( end.X() < 0 || end.X() > 256.35 || end.Y() > 116.5 || end.Y() < -116.5 || end.Z() < 0 || end.Z() > 1036.8 ){

	   std::cout<<"\nEVENT, "<<_event<<" Det energy: "<<det.E()<<std::endl ; //", "<<st.E()<<std::endl ;
	   //std::cout<<" ENd! " <<end.X()<<", "<<end.Y()<<", "<<end.Z()<<std::endl;
	 }
	 else{
	   shr_ids.emplace_back(si);
	 
	 }

	 
          //if ( s.DetProfile().E() > 0 ){
          //  shr_ids.emplace_back(si) ;
          // //if( det.X() < 0 || det.X() > 256.35 || det.Y() < -116.5 || det.Y() > 116.5 || det.Z() < 0 || det.Z() > 1036.8 ) 
	  // //  out = true; 
          //}
	  //else {
	  //  std::cout<<"Det energy: "<<det.E()<<", "<<st.E()<<std::endl ;
	  //  std::cout<<"det: "<<det.X()<<", "<<det.Y()<<", "<<det.Z()<<std::endl ;
	  //  //std::cout<<"st: "<<st.X()<<", "<<st.Y()<<", "<<st.Z()<<std::endl ;
	  //}
	}
     }

     //if ( out ) _out_of_vol++;
     if ( shr_ids.size() < 2 ) _thresh++;
     if ( shr_ids.size() == 3 ) _dalitz++ ;

     //if ( shr_ids.size() < 2 ) 
     //  _thresh++ ;
     //else if ( out ) 
     //  _out_of_vol++ ;
     //else if ( shr_ids.size() == 3 ) 
     //  _dalitz++ ;

    }

    return true;
  }

  bool AcceptanceStudy::finalize() {

   std::cout<<"Total true : "<<_infv_ccpi0<<std::endl ; 
   std::cout<<"Below Threshold: "<<_thresh<<std::endl ;
   std::cout<<"Dalitz:  "<<_dalitz<<std::endl ;
   std::cout<<"Shower out of volume :  "<<_out_of_vol<<std::endl ;
  
    return true;
  }

}
#endif
