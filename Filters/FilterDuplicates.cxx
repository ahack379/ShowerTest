#ifndef LARLITE_FILTERDUPLICATES_CXX
#define LARLITE_FILTERDUPLICATES_CXX

#include "FilterDuplicates.h"

namespace larlite {

  bool FilterDuplicates::initialize() {

    _map_v.resize(5,std::multimap<float,float>());

    _event = 0;
    _event_no_dup = 0;

    return true;
  }
  
  bool FilterDuplicates::analyze(storage_manager* storage) {
   
   _event++;

   auto _map = _map_v.at(storage->run_id());
   
   auto it = _map.find(storage->subrun_id());
   bool foundit = false;

   if( it != _map.end() ){
    while ( it->first == storage->subrun_id() ){  
      auto temp_event = it->second ; 
      if( temp_event == storage->event_id() )
        foundit = true;

      it++; 
      }   
     if ( !foundit)
      _map.emplace(storage->subrun_id(), storage->event_id() );
       
     else{
       std::cout<<"Duplicates "<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl;
       return false ;
     }
    }   
   else 
      _map.emplace(storage->subrun_id(), storage->event_id() );


    _event_no_dup++ ;

    return true;
  }

  bool FilterDuplicates::finalize() {

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
