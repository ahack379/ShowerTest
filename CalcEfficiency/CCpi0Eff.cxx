#ifndef LARLITE_CCPI0EFF_CXX
#define LARLITE_CCPI0EFF_CXX

#include "CCpi0Eff.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"

namespace larlite {

  bool CCpi0Eff::initialize() {    

    _event = 0; 
    _signal = 0;
    _event_list.clear();

    return true;
  }
  
  bool CCpi0Eff::analyze(storage_manager* storage) {

    //std::cout<<"\nEvent is : "<<_event <<", "<<storage->event_id()<<", "<<storage->subrun_id()<<std::endl ;

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


    _event++ ;

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
    auto e = traj.at(traj.size() - 1).E();

    if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
        return false;

    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0;
    int n_mu = 0;
    
    for ( auto const & p : parts ){
     
        if( p.StatusCode() == 1 && p.PdgCode() == 111 )
          n_pi0 ++;
        if( p.StatusCode() == 1 && p.PdgCode() == 13 )
          n_mu ++;
     }

    if ( n_pi0 == 1 && n_mu == 1 && e > 0.5) 
      _event_list.emplace_back(_event-1);

    return true;
  }

  bool CCpi0Eff::finalize() {

    std::cout<<"CCpi0 are "<<float(_event_list.size())/(_event)*100<<"\% of BNB ("<<_event_list.size()<<"/"<<_event<<")"<<std::endl ;

    std::cout<<"\n\n"<<_event_list.size()<<" in Event list :" <<std::endl ;
    for( auto const & e : _event_list) std::cout<<e<<", ";

    return true;
  }

}
#endif
