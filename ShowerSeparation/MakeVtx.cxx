#ifndef LARLITE_MAKEVTX_CXX
#define LARLITE_MAKEVTX_CXX

#include "MakeVtx.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool MakeVtx::initialize() {

    _id = -1;
    _event = 0;

    return true;
  }
  
  bool MakeVtx::analyze(storage_manager* storage) {

    _id++;  
     
    auto ev_vtx = storage->get_data<event_vertex>("mcvertex"); 

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
       else if( _mu ){
          auto ev_mct= storage->get_data<event_mctrack>("mcreco"); 
          if(!ev_mct || !ev_mct->size() ) return false;

          auto t = ev_mct->at(0);
          xyz[0] = t.Start().X() + _mu_offset ;
          xyz[1] = t.Start().Y();
          xyz[2] = t.Start().Z();
           }

	else if( _bnb ){
          auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
          if(!ev_mctruth || !ev_mctruth->size() ) return false;

          auto nu = ev_mctruth->at(0).GetNeutrino();
	  auto parts = ev_mctruth->at(0).GetParticles();

	  std::cout<<"\nEvent "<< _event<<", Interaction Type: "<<nu.InteractionType()<<std::endl; 
          //for ( auto const & p : parts ){
	  //  if( p.StatusCode() == 1 ){
	  //     std::cout<<"PDG : "<<p.PdgCode()<<std::endl ;
	  //     std::cout<<"Location : "<<p.Trajectory().at(0).X()<<", "<<p.Trajectory().at(0).Y()<<", "<<p.Trajectory().at(0).Z()<<std::endl ;
	  //     }
	  //  }

	  auto traj = nu.Nu().Trajectory();
          xyz[0] = traj.at(traj.size() - 1).X() + _bnb_offset ;
          xyz[1] = traj.at(traj.size() - 1).Y();
          xyz[2] = traj.at(traj.size() - 1).Z();

	  if( xyz[0] < 20 || xyz[0] > 236.5 || xyz[1] < -96.5 || xyz[1] > 96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
	    return false ;

	  //std::cout<<"And vertex location: "<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<std::endl;
	   }

        vertex new_vtx(xyz,_id) ;
        ev_vtx->push_back(new_vtx);

        }
	else return false;

	_event ++ ;


    return true;
  }

  bool MakeVtx::finalize() {

    return true;
  }

}
#endif
