#ifndef LARLITE_SEPARATECCNC_CXX
#define LARLITE_SEPARATECCNC_CXX

#include "SeparateCCNC.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

namespace larlite {

  bool SeparateCCNC::initialize() {    
    _signal = 0;
    _events = 0;

    return true;
  }
  
  bool SeparateCCNC::analyze(storage_manager* storage) {

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
      if(!ev_mctruth || !ev_mctruth->size() ) return false;

    auto ev_mct = storage->get_data<event_mctrack>("mcreco");

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    for ( auto const & t : *ev_mct ){
       if ( t.PdgCode() == 13 || t.PdgCode() == -13 ){

         auto st = t.Start() ;
         auto dist = sqrt( pow(xyz[0] - st.X(),2) + pow(xyz[1] - st.Y(),2) +pow(xyz[2] - st.Z(),2) );  

          if( dist < 0.4 ) { 

            if( _getNC ) return false;

            _signal++ ;
            //_event_list.emplace_back(_events);
            //std::cout<<"Foudn a CCpi0!!!!!!!!!!!"<<std::endl ;
            return true;
            }   
          }   
        }   

    if ( _getNC ){
     _signal++;
     return true;
     }
    else return false ;

  }

  bool SeparateCCNC::finalize() {

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
