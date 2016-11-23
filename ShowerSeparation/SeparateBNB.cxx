#ifndef LARLITE_SEPARATEBNB_CXX
#define LARLITE_SEPARATEBNB_CXX

#include "SeparateBNB.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/vertex.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"

namespace larlite {

  bool SeparateBNB::initialize() {    

    _event = 0; 
    return true;
  }
  
  bool SeparateBNB::analyze(storage_manager* storage) {

    std::cout<<"\nEvent is : "<<_event <<std::endl ;
    _event++ ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) return false;

    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    auto ev_mct = storage->get_data<event_mctrack>("mcreco");

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();

    // Some constants 
    double det_drift_velocity = ::larutil::LArProperties::GetME()->DriftVelocity(); ///< cm/us
    double event_time = ( traj.at(traj.size()-1).T() - 3200 ) * 0.5 ; // ticks * us / tick = us
    double shift_x = event_time * det_drift_velocity ; //cm

    xyz[0] = traj.at(traj.size() - 1).X(); // + shift_x;
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    auto const& geomH = ::larutil::GeometryHelper::GetME();
    auto vtxWT  = geomH->Point_3Dto2D(xyz,2);
    
    auto vtx_w = vtxWT.w / 0.3; //geomH->WireToCm();
    auto vtx_t = vtxWT.t / 0.05 + 800; // geomH->TimeToCm();

    //std::cout<<"X Y Z at interaction is  : "<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<std::endl ;
    //std::cout<<nu.CCNC()<<", Wire and time for vertex : "<<vtx_w<<", "<<vtx_t<<std::endl ;
    //std::cout << "Run   " << storage->run_id() << std::endl;
    //std::cout << "Subrun " << storage->subrun_id() << std::endl;
    //std::cout << "Event " << storage->event_id() << std::endl;

    bool pi0 = false ;

    //auto const & t = nu.InteractionType();

    //if( (t == 1004 || t == 1011 || t == 1080 || t == 1086 || t == 1090) )  
    //  std::cout<<"CC pi0 actual "<<std::endl ;

    //for ( auto const & s : *ev_mcs ){

    //   if ( s.MotherPdgCode() == 111 ){ 
    //     auto st = s.Start() ;
    //     auto dist = sqrt( pow(xyz[0] - st.X(),2) + pow(xyz[1] - st.Y(),2) +pow(xyz[2] - st.Z(),2) ); 
    //     double xyz2[3] = {0.};
    //     xyz2[0] = st.X(); 
    //     xyz2[1] = st.Y();
    //     xyz2[2] = st.Z();

    //     auto vtx2  = geomH->Point_3Dto2D(xyz2,2);

    //     //std::cout<<"X Y Z at interaction is  : "<<xyz2[0]<<", "<<xyz2[1]<<", "<<xyz2[2]<<", E: "
    //     //        <<s.DetProfile().E()<<", Mother: "<<s.MotherPdgCode()<<std::endl ;
    //     //std::cout<<"PI0 Wire and time: "<<vtx2.w/0.3<<", "<<vtx2.t/0.05+800<<", Mom: "<<s.MotherPdgCode()<<std::endl ;

    //      if( dist < 0.01 ) {
    //        std::cout<<"FOUND!"<<dist <<std::endl;
    //        pi0 = true;
    //        break;
    //       }    
    //     }     
    //   }     

    //if ( pi0 ){
    //  if( _get_pi0s ) return true;
    //  else return false;
    // }
    //else{
    //  if( _get_pi0s ) return false ;
    //  else return true ;
    //  }

      auto parts = ev_mctruth->at(0).GetParticles();
      
      for ( auto const & p : parts ){
         //if (p.PdgCode() == 111 )
         //  std::cout<<"status "<<p.StatusCode()<<std::endl ;
        
        if( p.StatusCode() == 1 && p.PdgCode() == 111 ) {
          std::cout<<"FOUND PI0!! "<<std::endl;
	  pi0 = true;
	  break;
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

  }

  bool SeparateBNB::finalize() {
  
    return true;
  }

}
#endif
