#ifndef LARLITE_EVTWEIGHTFILTER_CXX
#define LARLITE_EVTWEIGHTFILTER_CXX

#include "EvtWeightFilter.h"
#include <sstream>

namespace larlite {

  bool EvtWeightFilter::initialize() {

    _event = 0;
    _pi0_map.resize(10,std::multimap<float,float>());
    _map_v.resize(10,std::multimap<float,float>());

    // load in all run subrun and event ids from the two or one shower selections (run over full mcbnbcos) 
    _file.open("vtx_info.txt",std::ios_base::in);

    for(std::string line; std::getline(_file, line); )   //read stream line by line
    {
        std::istringstream in(line);      //make a stream for the line itself

        float run;
        in >> run;                  //and read the first whitespace-separated token
        float subrun;
        in >> subrun;                  //and read the first whitespace-separated token
	    float event ;
	    in >> event ;

	    _pi0_map.at(run).emplace(subrun,event);
    }

    return true;
  }
  
  // Here we're looping over events in the EVENTWEIGHT output trying to filter out only the events we have also in the mcbnbcos sample
  bool EvtWeightFilter::analyze(storage_manager* storage) {

    //std::cout<<"GO Subrun and event : "<<storage->subrun_id()<<", "<<storage->event_id() <<std::endl ;
    _event ++;

    // Check 2 things about this event in this analyze funtion:
    // 1) First, verify that the current eventweight event made it to the final two (or one) shower sample 
    auto r = storage->run_id() ;
    bool foundit = false;
    auto it = _pi0_map.at(r).find(storage->subrun_id());

    if( it != _pi0_map.at(r).end() ){
      while ( it->first == storage->subrun_id() ){  
        auto temp_event = it->second ; 
        if( temp_event == storage->event_id() ){
          foundit = true;
          break;
        }   
        it++; 
      }
    }

    if ( !foundit ) return false ;

    // 2) Second, verify that the current eventweight event is not a duplicate of an event already identified in number 1).
    foundit = false;
    auto it2 = _map_v.at(r).find(storage->subrun_id());

    if( it2 != _map_v.at(r).end() ){
     while ( it2->first == storage->subrun_id() ){   
       auto temp_event = it2->second ; 
       if( temp_event == storage->event_id() )
         foundit = true;

       it2++; 
       }    
      if ( !foundit)
       _map_v.at(r).emplace(storage->subrun_id(), storage->event_id() );
      else{
        std::cout<<"Event: "<<_event<<std::endl;
        std::cout<<"Duplicates "<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl;
        return false ;
      }    
    }    
    else 
      _map_v.at(r).emplace(storage->subrun_id(), storage->event_id() );



    std::cout<<"Found one! "<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl; 

    return true;
  }

  bool EvtWeightFilter::finalize() {

    return true;
  }

}
#endif
