#ifndef LARLITE_NOMES_SEPARATEBNB_CXX
#define LARLITE_NOMES_SEPARATEBNB_CXX

#include "NoMes_SeparateBNB.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/vertex.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"

namespace larlite {

  bool NoMes_SeparateBNB::initialize() {    

    _event = 0; 
    _signal = 0;
    return true;
  }
  
  bool NoMes_SeparateBNB::analyze(storage_manager* storage) {

    std::cout<<"\nEvent is : "<<_event <<std::endl ;
    _event++ ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) return false;

    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();

    // Some constants 
    double det_drift_velocity = ::larutil::LArProperties::GetME()->DriftVelocity(); ///< cm/us
    double event_time = ( traj.at(traj.size()-1).T() - 3200 ) * 0.5 ; // ticks * us / tick = us

    xyz[0] = traj.at(traj.size() - 1).X(); // + shift_x;
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    auto const& geomH = ::larutil::GeometryHelper::GetME();
    auto vtxWT  = geomH->Point_3Dto2D(xyz,2);
    
    //std::cout<<"X Y Z at interaction is  : "<<xyz[0]<<", "<<xyz[1]<<", "<<xyz[2]<<std::endl ;
    //std::cout<<nu.CCNC()<<", Wire and time for vertex : "<<vtx_w<<", "<<vtx_t<<std::endl ;
    //std::cout << "Run   " << storage->run_id() << std::endl;
    //std::cout << "Subrun " << storage->subrun_id() << std::endl;
    //std::cout << "Event " << storage->event_id() << std::endl;

    bool pi0 = false ;

    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0;
    int n_mes = 0;
    int n_lep = 0;
    int n_mu = 0;
      
    for ( auto const & p : parts ){

     if( p.StatusCode() == 1 && p.PdgCode() == 111 ){
        n_pi0 += 1;
        }   

      if( p.StatusCode() == 1 && p.PdgCode() == 13 )
        n_mu += 1;

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 211 || abs(p.PdgCode()) == 321 ||  p.PdgCode() == 130 
                                   || p.PdgCode() == 310 || abs(p.PdgCode()) == 311 ) ){
        n_mes += 1;
        break ;
        }   

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 11 || p.PdgCode() == 22 )){ 
        n_lep += 1;
        break ;
        }   
      }   

      if( n_mu == 1 && n_pi0 == 1 && n_lep == 0 && n_mes == 0){ 
        pi0 = true;
        _signal++;
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

  bool NoMes_SeparateBNB::finalize() {
  
    std::cout<<"OUT SIG! "<<_signal<<std::endl ;

    return true;
  }

}
#endif
