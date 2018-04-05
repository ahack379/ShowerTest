#ifndef LARLITE_WRITEWEIGHTS_CXX
#define LARLITE_WRITEWEIGHTS_CXX

#include "WriteWeights.h"
#include "DataFormat/vertex.h"
#include "DataFormat/track.h"
#include "DataFormat/mceventweight.h"

namespace larlite {

  bool WriteWeights::initialize() {

    //event_file.open("GenieWeightStudy.txt");
    _event = 0;

    return true;
  }
  
  bool WriteWeights::analyze(storage_manager* storage) {
  

    //std::cout<<"*****************************************Event: "<<_event <<std::endl; 
    _event ++ ;
    // Run this module on the final sample of pi0 info.  
    auto run = storage->run_id();
    auto event = storage->event_id();
    auto subrun = storage->subrun_id();

    auto ev_wgt = storage->get_data<event_mceventweight>(_event_producer); 
    //auto ev_wgt_flux = storage->get_data<event_mceventweight>("fluxeventweight");

    if( !ev_wgt || !ev_wgt->size() ){
      std::cout<<"No event weights..." <<std::endl;
      return false;
    }

    //event_file << run << " " << subrun << " " << event <<" "; //<< " "
    std::cout<<run<<" "<<subrun<<" "<<event<<" " ; //<<std::endl ;

    auto wgt  = ev_wgt->at(0).GetWeights();
    auto w_v = wgt.begin()->second; //

    //std::cout<<"Number of weights : "<<ev_wgt->size()<<std::endl ;

    std::vector<double> weights_v;

    for ( auto const & m : wgt ) { 
       if (m.first == "bnbcorrection_FluxHist" ) {
           continue;
       }
       //std::cout<<"Parameter: "<<m.first<<", "<<m.second.size() <<std::endl;

       for ( auto const & w : m.second){
         //event_file<<w<<" " ;
         std::cout<<w<<" "; //<<std::endl ;
        }   
    }   

    std::cout<<"\n";
    //event_file<<"\n";


    return true;
  }

  bool WriteWeights::finalize() {

    return true;
  }

}
#endif
