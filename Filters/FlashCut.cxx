#ifndef LARLITE_FLASHCUT_CXX
#define LARLITE_FLASHCUT_CXX

#include "FlashCut.h"
#include "DataFormat/opflash.h"

namespace larlite {

  bool FlashCut::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool FlashCut::analyze(storage_manager* storage) {

    auto ev_flash = storage->get_data<event_opflash>("opflashSat");

    if ( !ev_flash || ev_flash->size() == 0 ) return false ;

    for( auto const & f : *ev_flash ){
      if( f.TotalPE() > 50 && f.Time() >3.65 && f.Time() < 5.25)
        return true;
        //std::cout<<"Test. "<<std::endl;
       }
  
    return false;
  }

  bool FlashCut::finalize() {

    return true;
  }

}
#endif
