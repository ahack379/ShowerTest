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
        //std::cout<<"Duplicates "<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl;
        return false ;
      }     
     }    
    else 
       _map_v.at(r).emplace(storage->subrun_id(), storage->event_id() );

    std::cout<<storage->run_id()<<" "<<storage->subrun_id()<<" "<<storage->event_id()<<std::endl ;

    evt++; 
  
    return true;
  }

  bool GetEventRunNumbers::finalize() {

   //for ( auto & e : _map )
   // std::cout<<e.first<<" "<<e.second<<std::endl ;

    return true;
  }

}
#endif
