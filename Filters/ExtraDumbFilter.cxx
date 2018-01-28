#ifndef LARLITE_EXTRADUMBFILTER_CXX
#define LARLITE_EXTRADUMBFILTER_CXX

#include "ExtraDumbFilter.h"

namespace larlite {

  void ExtraDumbFilter::AddSubrunEvent(int subrun, int event) {
    _evt_m[subrun] = event ;
  }

  bool ExtraDumbFilter::initialize() {


    return true;
  }
  
  bool ExtraDumbFilter::analyze(storage_manager* storage) {
 
    if ( _evt_m.find(storage->subrun_id()) != _evt_m.end() ){
      if ( _evt_m[int(storage->subrun_id())] == int(storage->event_id()) )
        return true;
      else return false;
    }
    //if ( storage->subrun_id() == _subrun && storage->event_id() == _event)
    //  return true;
  
    return false;
  }

  bool ExtraDumbFilter::finalize() {
  
    return true;
  }

}
#endif
