#ifndef LARLITE_GETEVENTRUNNUMBERS_CXX
#define LARLITE_GETEVENTRUNNUMBERS_CXX

#include "GetEventRunNumbers.h"
#include <sstream>

namespace larlite {

  bool GetEventRunNumbers::initialize() {
    
    evt = 0;

    return true;
  }
  
  bool GetEventRunNumbers::analyze(storage_manager* storage) {


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

    }
   else 
      _map.emplace(storage->subrun_id(), storage->event_id() );

   // //std::cout<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl ;

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
