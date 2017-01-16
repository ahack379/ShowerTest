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

    //std::cout<<"\nEvent is : "<<_event <<std::endl ;
    _event++ ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
      }

    //auto ev_mctruth_cos = storage->get_data<event_mctruth>("corsika"); 
    //if(!ev_mctruth_cos || !ev_mctruth_cos->size() ) {
    //  std::cout<<"Event has no mctruth corsika info "<<std::endl;
    //  return false;
    //  }

    auto ev_mcshr= storage->get_data<event_mcshower>("mcreco"); 
    if(!ev_mcshr || !ev_mcshr->size() ) {
      std::cout<<"Event has no mcshower info "<<std::endl;
      return false;
      }

    //for ( auto const & s : *ev_mcshr ){
    //  if ( abs(s.MotherPdgCode()) != 13 && s.DetProfile().E() > 0 && s.Origin() == 1)
    //    std::cout<<"PDG: "<<s.MotherPdgCode()<<", "<<s.Start().E()<<std::endl ;
    //  }


    auto parts = ev_mctruth->at(0).GetParticles();
    bool pi0 = false;
    bool mu  = false ;

    //std::cout<<"N particles: "<<parts.size()<<std::endl ;
    float energy = 0;

    for ( auto const & p : parts ){

      //std::cout<<"Each pdg: "<<p.PdgCode()<<", "<<p.Process()<<std::endl ;
      
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){ 
        pi0 = true;
        energy = p.Trajectory().at(0).E();
	//std::cout<<"Found pi0 "<<std::endl ;
	}

      if( p.StatusCode() == 1 && p.PdgCode() == 13 ){
        mu = true; 
	//std::cout<<"Found mu "<<std::endl ;
	}
        
      if( mu && pi0 ){
        std::cout<<"**********************FOUND CCPI0!! "<<energy*1000<<std::endl;
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
    for( auto const & e : _event_list ) std::cout<<e<<", ";


  
    return true;
  }

}
#endif
