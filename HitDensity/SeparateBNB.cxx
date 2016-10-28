#ifndef LARLITE_SEPARATEBNB_CXX
#define LARLITE_SEPARATEBNB_CXX

#include "SeparateBNB.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool SeparateBNB::initialize() {    

    return true;
  }
  
  bool SeparateBNB::analyze(storage_manager* storage) {

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
      if(!ev_mctruth || !ev_mctruth->size() ) return false;

      auto parts = ev_mctruth->at(0).GetParticles();
      
      for ( auto const & p : parts ){
        
        if( p.PdgCode() == 11 || p.PdgCode() == 22 || p.PdgCode() == 111 ){
          if( _get_shower_events )
            return true ;
          else
            return false ;
           }
         }
     
    if( _get_shower_events )
      return false;
    else
      return true;
  }

  bool SeparateBNB::finalize() {

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
