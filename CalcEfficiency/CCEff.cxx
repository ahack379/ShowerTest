#ifndef LARLITE_CCEFF_CXX
#define LARLITE_CCEFF_CXX

#include "CCEff.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"

namespace larlite {

  bool CCEff::initialize() {    

    _event = 0; 
    _signal = 0;

    _n_cc = 0;    // 0 

    return true;
  }
  
  bool CCEff::analyze(storage_manager* storage) {

    _event++ ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
      }

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();
    auto e = traj.at(traj.size() - 1).E();

    bool infv = true;

    if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
        infv = false;

    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0;
    int n_mu = 0;
    
    //for ( auto const & p : parts ){
    // 
    //    if( p.StatusCode() == 1 && p.PdgCode() == 13 )
    //      n_mu ++;
    // }

     if( nu.CCNC() == 0 && nu.Nu().PdgCode() == 14 && e > 0.5 && infv){ 
       _event_list.emplace_back(_event-1);
       _n_cc++; 
     }


    return true;
  }

  bool CCEff::finalize() {

    std::cout<<"Total accounted backgrounds: "<< _n_cc<<std::endl;

    return true;
  }

}
#endif
