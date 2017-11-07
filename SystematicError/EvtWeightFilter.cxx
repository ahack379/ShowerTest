#ifndef LARLITE_EVTWEIGHTFILTER_CXX
#define LARLITE_EVTWEIGHTFILTER_CXX

#include "EvtWeightFilter.h"
#include <sstream>

namespace larlite {

  bool EvtWeightFilter::initialize() {

    _file.open("vtx_info.txt",std::ios_base::in);
    //_file.open("mcc8_pi0cuts.txt",std::ios_base::in);
    //_file.open("eventweightCheck.txt",std::ios_base::in);

    for(std::string line; std::getline(_file, line); )   //read stream line by line
    {
        std::istringstream in(line);      //make a stream for the line itself
    
        float subrun;
        in >> subrun;                  //and read the first whitespace-separated token
	float event ;
	in >> event ;
	_pi0_map.emplace(subrun,event);

    }

    return true;
  }
  
  bool EvtWeightFilter::analyze(storage_manager* storage) {

    //std::cout<<"GO Subrun and event : "<<storage->subrun_id()<<", "<<storage->event_id() <<std::endl ;

    bool foundit = false;
    auto it = _pi0_map.find(storage->subrun_id());

    if( it != _pi0_map.end() ){
      while ( it->first == storage->subrun_id() ){ 
        auto temp_event = it->second ; 
        if( temp_event == storage->event_id() ){
          foundit = true;
	  break ;
	  }

        it++; 
         }
    }
   if ( !foundit) return false;

    return true;
  }

  bool EvtWeightFilter::finalize() {

    return true;
  }

}
#endif
