#ifndef LARLITE_RUNSUBRUNEVENT_CXX
#define LARLITE_RUNSUBRUNEVENT_CXX

#include "RunSubrunEvent.h"
#include "DataFormat/vertex.h"

namespace larlite {

  bool RunSubrunEvent::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool RunSubrunEvent::analyze(storage_manager* storage) {

    //auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");
    //auto vtx = ev_vtx->at(0);
    
    std::cout<<storage->run_id()<<" "<<storage->subrun_id()<<" "<<storage->event_id()<<std::endl ; //" "<<vtx.X()<<" "<<vtx.Y()<<" "<<vtx.Z()<<std::endl ;
    //std::cout<<storage->subrun_id()<<" "<<storage->event_id()<<std::endl ;
  
    return true;
  }

  bool RunSubrunEvent::finalize() {

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
    // if(_fout) { _fout->cd(); h1->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    return true;
  }

}
#endif
