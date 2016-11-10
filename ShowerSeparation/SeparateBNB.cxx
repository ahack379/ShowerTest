#ifndef LARLITE_SEPARATEBNB_CXX
#define LARLITE_SEPARATEBNB_CXX

#include "SeparateBNB.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

namespace larlite {

  bool SeparateBNB::initialize() {    

    return true;
  }
  
  bool SeparateBNB::analyze(storage_manager* storage) {

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
      if(!ev_mctruth || !ev_mctruth->size() ) return false;

    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    auto ev_mct = storage->get_data<event_mctrack>("mcreco");

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    bool pi0 = false ;

    for ( auto const & s : *ev_mcs ){

       if ( s.MotherPdgCode() == 111 ){

         auto st = s.Start() ;
         auto dist = sqrt( pow(xyz[0] - st.X(),2) + pow(xyz[1] - st.Y(),2) +pow(xyz[2] - st.Z(),2) ); 
          if( dist < 0.4 ) {
            pi0 = true;
            break;
           }    
         }     
       }     

    if ( pi0 ){
      if( _get_pi0s ) return true;
      else return false;
     }
    else{
      if( _get_pi0s ) return false ;
      else return true ;
      }

//      auto parts = ev_mctruth->at(0).GetParticles();
//      
//      for ( auto const & p : parts ){
//        
//        if( p.StatusCode() == 1 && ( p.PdgCode() == 111 ) ){
//          if( _get_shower_events )
//            return true ;
//          else
//            return false ;
//           }
//         }
//     
//    if( _get_shower_events )
//      return false;
//    else
//      return true;
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
