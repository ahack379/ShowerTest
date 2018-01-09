#ifndef LARLITE_EXTRADUMBFILTER_CXX
#define LARLITE_EXTRADUMBFILTER_CXX

#include "ExtraDumbFilter.h"

namespace larlite {

  bool ExtraDumbFilter::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool ExtraDumbFilter::analyze(storage_manager* storage) {
 
    if ( storage->subrun_id() == _subrun && storage->event_id() == _event){
      return true;
    }
  
    return false;
  }

  bool ExtraDumbFilter::finalize() {
  
    return true;
  }

}
#endif
