#ifndef LARLITE_DATAGETEVENTRUNNUMBERS_CXX
#define LARLITE_DATAGETEVENTRUNNUMBERS_CXX

#include "DataGetEventRunNumbers.h"
#include <sstream>

namespace larlite {

  bool DataGetEventRunNumbers::initialize() {
    
    evt = 0;

    return true;
  }
  
  bool DataGetEventRunNumbers::analyze(storage_manager* storage) {


   auto it = _map.find(storage->subrun_id());
   bool foundit = false;

   auto tmp = std::make_pair(storage->subrun_id(),storage->event_id());

   //if( it != _map.end() ){
   // while ( it->first == storage->subrun_id() ){ 
   //   auto temp_event = it->second ; 
   //   if( temp_event == storage->event_id() )
   //     foundit = true;

   //   it++; 
   //   }
   //  if ( !foundit)
   //   _map.emplace(storage->run_id(), tmp ); 

   // }
   //else 
   //   _map.emplace(storage->run_id(),tmp) ; 

    //std::cout<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl ;

   _map.emplace(storage->run_id(),tmp) ; 
    //_map.emplace(storage->subrun_id(), storage->event_id() );

    evt++; 
  
    return true;
  }

  bool DataGetEventRunNumbers::finalize() {

   for ( auto & e : _map )
    std::cout<<e.first<<" "<<(e.second).first<<" "<<(e.second).second<<std::endl ;

    return true;
  }

}
#endif
