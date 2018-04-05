#ifndef LARLITE_GETEVENTRUNNUMBERS_CXX
#define LARLITE_GETEVENTRUNNUMBERS_CXX

#include "GetEventRunNumbers.h"
#include "DataFormat/mctruth.h"
#include <sstream>

namespace larlite {

  bool GetEventRunNumbers::initialize() {
    
    evt = 0;
    _map_v.resize(10,std::multimap<float,float>());

    return true;
  }
  
  bool GetEventRunNumbers::analyze(storage_manager* storage) {

    auto r = storage->run_id() ;
    auto it = _map_v.at(r).find(storage->subrun_id());
    bool foundit = false;

    if( it != _map_v.at(r).end() ){
     while ( it->first == storage->subrun_id() ){  
       auto temp_event = it->second ; 
       if( temp_event == storage->event_id() )
         foundit = true;

       it++; 
       }   
      if ( !foundit)
       _map_v.at(r).emplace(storage->subrun_id(), storage->event_id() );
    
      else{
        //std::cout<<"Event: "<<_event<<std::endl;
        std::cout<<"Duplicates "<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl;
        return false ;
      }     
     }    
    else 
       _map_v.at(r).emplace(storage->subrun_id(), storage->event_id() );

    std::cout<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl ;

    //auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    //if(!ev_mctruth || !ev_mctruth->size() ) { 
    //  std::cout<<"Event has no mctruth info "<<std::endl;
    //  return false;
    //  }   

    //auto & truth = ev_mctruth->at(0);
    //auto & nu  = truth.GetNeutrino();

    //double xyz[3] = {0.};
    //auto traj = nu.Nu().Trajectory();
    //xyz[0] = traj.at(traj.size() - 1).X();
    //xyz[1] = traj.at(traj.size() - 1).Y();
    //xyz[2] = traj.at(traj.size() - 1).Z();
    //auto e = traj.at(traj.size() - 1).E();

    //bool infv = true;

    //if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
    //    infv = false ;

    //auto parts = ev_mctruth->at(0).GetParticles();
    //int n_pi0 = 0;
    //int n_mu = 0;
    //
    //for ( auto const & p : parts ){
    //
    //    if( p.StatusCode() == 1 && p.PdgCode() == 111 )
    //      n_pi0 ++; 
    //    if( p.StatusCode() == 1 && p.PdgCode() == 13 )
    //      n_mu ++; 
    // }   

    //if ( n_pi0 == 1 && n_mu == 1 && infv){
    //  _map.emplace(storage->subrun_id(),storage->event_id());

    //} 
    //  _event_list.emplace_back(_event-1);



   //_map.emplace(storage->subrun_id(), storage->event_id() );

    evt++; 
  
    return true;
  }

  bool GetEventRunNumbers::finalize() {

   for ( auto & e : _map )
    std::cout<<e.first<<" "<<e.second<<std::endl ;

    return true;
  }

}
#endif
