#ifndef LARLITE_SEL2CCPI0EFF_CXX
#define LARLITE_SEL2CCPI0EFF_CXX

#include "Sel2CCpi0Eff.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool Sel2CCpi0Eff::initialize() {

  std::cout<<"Entering initialize "<<std::endl ;

    _events = 0 ;
    _signal = 0;

    return true;
  }
  
  bool Sel2CCpi0Eff::analyze(storage_manager* storage) {

    
    // This is a calculation of CCpi0 in the BNB file 
    // Can also be run on sel2 filter to calculate the
    // eff of sel2 filter for OUR signal
    auto ev_truth = storage->get_data<event_mctruth>("generator");
    auto & truth = ev_truth->at(0);
    auto & nu  = truth.GetNeutrino();
    
    auto t = nu.InteractionType();

    //These are the Interaction codes taht inclides pi0 and muon
    if( t == 1004 || t == 1011 || t == 1080 || t == 1086 || t == 1090 ) 
      _signal++ ;

    _events++;

    return true;
  }

  bool Sel2CCpi0Eff::finalize() {

    std::cout<<"CCpi0 are "<<float(_signal)/_events*100<<"\% of BNB ("<<_signal<<"/"<<_events<<")"<<std::endl ;
  
    return true;
  }

}
#endif
