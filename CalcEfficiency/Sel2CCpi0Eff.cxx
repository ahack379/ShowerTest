#ifndef LARLITE_SEL2CCPI0EFF_CXX
#define LARLITE_SEL2CCPI0EFF_CXX

#include "Sel2CCpi0Eff.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool Sel2CCpi0Eff::initialize() {

  std::cout<<"Entering initialize "<<std::endl ;

    _events = 0 ;
    _signal = 0;

    _event_list.clear();

    return true;
  }
  
  bool Sel2CCpi0Eff::analyze(storage_manager* storage) {

    
    // This is a calculation of CCpi0 in the BNB file 
    // Can also be run on sel2 filter to calculate the
    // eff of sel2 filter for OUR signal
    auto ev_truth = storage->get_data<event_mctruth>("generator");
    if(!ev_truth || !ev_truth->size() ) return false;

    auto & truth = ev_truth->at(0);
    auto & nu  = truth.GetNeutrino();
    //auto t = nu.InteractionType();

    //These are the Interaction codes that include pi0 and muon
    //if( t == 1004 || t == 1011 || t == 1080 || t == 1086 || t == 1090 ){ 

    double xyz[3] = {0.};

    auto parts = ev_truth->at(0).GetParticles();

    for ( auto const & p : parts ){

      if( p.StatusCode() == 1 && ( p.PdgCode() == 111 ) ){

          auto traj = nu.Nu().Trajectory();
          xyz[0] = traj.at(traj.size() - 1).X();
          xyz[1] = traj.at(traj.size() - 1).Y();
          xyz[2] = traj.at(traj.size() - 1).Z();

          auto ptraj = p.Trajectory() ;
          if( !ptraj.size() ) continue;

          auto dist = sqrt( pow(xyz[0] - ptraj.at(0).X(),2) + pow(xyz[1] - ptraj.at(0).Y(),2) +pow(xyz[2] - ptraj.at(0).Z(),2) ); 
         
          if( dist < 0.4 ) {
            _signal++ ;
            _event_list.emplace_back(_events);
            continue;
           }
        }

      }

    _events++;

    return true;
  }

  bool Sel2CCpi0Eff::finalize() {

    std::cout<<"CCpi0 are "<<float(_signal)/_events*100<<"\% of BNB ("<<_signal<<"/"<<_events<<")"<<std::endl ;

    std::cout<<"Event list :" <<std::endl ;
    for( auto const & e : _event_list ){
      
      std::cout<<e<<", ";
    
      }
  
    return true;
  }

}
#endif
