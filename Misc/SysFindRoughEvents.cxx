#ifndef LARLITE_SYSFINDROUGHEVENTS_CXX
#define LARLITE_SYSFINDROUGHEVENTS_CXX

#include "SysFindRoughEvents.h"

namespace larlite {

  bool SysFindRoughEvents::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    _file.open("sel2_ext.txt",std::ios_base::in);

    for(std::string line; std::getline(_file, line); )   //read stream line by line
    {   
        std::istringstream in(line);      //make a stream for the line itself
    
        float subrun;
        in >> subrun;                  //and read the first whitespace-separated token
        float event ;
        in >> event ;

        _sel2_ext_m[subrun] = event;
    }   
   
    _file.close();

    _file.open("pi0_ext.txt",std::ios_base::in);

    for(std::string line; std::getline(_file, line); )   //read stream line by line
    {   
        std::istringstream in(line);      //make a stream for the line itself
    
        float subrun;
        in >> subrun;                  //and read the first whitespace-separated token
        float event ;
        in >> event ;
        _pi0_ext_m[subrun] = event;
    }   
   
    _file.close();

    return true;
  }
  
  bool SysFindRoughEvents::analyze(storage_manager* storage) {
  
    auto subrun = storage->subrun_id();
    auto event = storage->event_id();

    if( _pi0_ext_m.find(subrun) == _pi0_ext_m.end() ){
      // At this point, we know the event passed cv pi0 cuts but failed enhancedext pi0 cuts
      if( _sel2_ext_m.find(subrun) != _sel2_ext_m.end() ){
        // At this point, we know the event also passed ext sel2 cuts
        if ( _sel2_ext_m[subrun] == event ){
          // Found the events we're interested in!
          std::cout<<subrun<<", "<<event<<std::endl ;
        }       
      }
    }
    // If this doesn't work, also consider events with matching subruns, but different event numbers
  
    return true;
  }

  bool SysFindRoughEvents::finalize() {

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
