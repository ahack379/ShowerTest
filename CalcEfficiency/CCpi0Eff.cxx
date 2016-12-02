#ifndef LARLITE_CCPI0EFF_CXX
#define LARLITE_CCPI0EFF_CXX

#include "CCpi0Eff.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
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

    std::cout<<"\nEvent is : "<<_event <<std::endl ;
    _event++ ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) return false;

    auto parts = ev_mctruth->at(0).GetParticles();
    bool pi0 = false;
    bool mu  = false ;
      
    for ( auto const & p : parts ){
      
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){ 
        pi0 = true;
	//std::cout<<"Found pi0 "<<std::endl ;
	}

      if( p.StatusCode() == 1 && abs(p.PdgCode()) == 13 ){
        mu = true; 
	//std::cout<<"Found mu "<<std::endl ;
	}
        
      if( mu && pi0 ){
        std::cout<<"**********************FOUND CCPI0!! "<<std::endl;
        _signal++;
	_event_list.emplace_back(_event-1);
        break;
        }
      }

    return true;
  }

  bool CCpi0Eff::finalize() {

     std::cout<<"CCpi0 are "<<float(_signal)/(_event)*100<<"\% of BNB ("<<_signal<<"/"<<_event<<")"<<std::endl ;

     std::cout<<_event_list.size()<<" in Event list :" <<std::endl ;
    //for( auto const & e : _event_list ) std::cout<<e<<", ";


  
    return true;
  }

}
#endif
