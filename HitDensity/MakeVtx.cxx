#ifndef LARLITE_MAKEVTX_CXX
#define LARLITE_MAKEVTX_CXX

#include "MakeVtx.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

namespace larlite {

  bool MakeVtx::initialize() {

    _id = -1;

    return true;
  }
  
  bool MakeVtx::analyze(storage_manager* storage) {

    _id++;  
     
    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex"); 

    storage->set_id(storage->run_id(), storage->subrun_id(), storage->event_id());

    if(ev_vtx){

      ev_vtx->reserve(1);
      double xyz[3] = {0.};

      if( _pi0 ){
        auto ev_mcs= storage->get_data<event_mcshower>("mcreco"); 
        if(!ev_mcs || !ev_mcs->size() ){std::cout<<"NO SHOW "<<std::endl ; return false;}
        auto s = ev_mcs->at(0);
          
         xyz[0] = s.Start().Position().X()+ _pi0_offset;
         xyz[1] = s.Start().Position().Y();
         xyz[2] = s.Start().Position().Z();
           }
       else{
          auto ev_mct= storage->get_data<event_mctrack>("mcreco"); 
          if(!ev_mct || !ev_mct->size() ) return false;

          auto t = ev_mct->at(0);
          xyz[0] = t.Start().X() + _mu_offset ;
          xyz[1] = t.Start().Y();
          xyz[2] = t.Start().Z();
           }

        vertex new_vtx(xyz,_id) ;
        ev_vtx->push_back(new_vtx);

        }
	else return false;


    return true;
  }

  bool MakeVtx::finalize() {

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
