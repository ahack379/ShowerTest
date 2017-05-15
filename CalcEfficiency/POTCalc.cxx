#ifndef LARLITE_POTCALC_CXX
#define LARLITE_POTCALC_CXX

#include "POTCalc.h"
#include "DataFormat/potsummary.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"

namespace larlite {

  bool POTCalc::initialize() {    

    _event = 0; 
    _event_list.clear();

    _tot_pot = 0;
    
     return true;
  }
  
  bool POTCalc::analyze(storage_manager* storage) {

   //auto it = _map.find(storage->subrun_id());
   //bool foundit = false;

   // if( it != _map.end() ){
   //   while ( it->first == storage->subrun_id() ){  
   //     auto temp_event = it->second ; 
   //     if( temp_event == storage->event_id() )
   //       foundit = true;

   //     it++; 
   //     }   
   //    if ( !foundit )
   //     _map.emplace(storage->subrun_id(), storage->event_id() );
   //    else return false;
   //   }   
   //  else 
   //     _map.emplace(storage->subrun_id(), storage->event_id() );

    _event++ ;

    // Now calculate the total POT + total numu neutrinos 
    auto ev_pot = storage->get_subrundata<potsummary>("generator"); 

    if( storage->subrun_id() != storage->last_subrun_id() )
      _tot_pot += ev_pot->totgoodpot ;

    return true;
  }

  bool POTCalc::finalize() {

    std::cout<<"Total POT: "<<_tot_pot<<std::endl;
    std::cout<<"Total Event: "<<_event<<std::endl;
  
    return true;
  }

}
#endif
