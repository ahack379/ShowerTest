#ifndef LARLITE_CROSSSECTION_CXX
#define LARLITE_CROSSSECTION_CXX

#include "CrossSection.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/vertex.h"
#include "DataFormat/potsummary.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mcflux.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"

namespace larlite {

  bool CrossSection::initialize() {    

    _event = 0; 
    _event_list.clear();
    
    _signal = 0;
    _tot_pot = 0.;
    _tot_event_in_AV = 0. ;
    _n_numu = 0; 
    _n_nu_all = 0;

    _xmin = 20. ;
    _xmax = 236.35 ;
    _ymin = -96.5 ;
    _ymax = 96.5 ;
    _zmin = 10. ;
    _zmax = 1026.8 ;
    
    _mean_e = 0.;

    //_xmin = 0. ;
    //_xmax = 256.35. ;
    //_ymin = -116.5 ;
    //_ymax = 116.5 ;
    //_zmin = 0. ;
    //_zmax = 1036.35 ;

    return true;
  }
  
  bool CrossSection::analyze(storage_manager* storage) {

    auto ev_temp = storage->event_id();
    auto subrun_temp = storage->subrun_id();

    std::cout<<"\nEvent is : "<<_event<<", "<<ev_temp<<", "<<subrun_temp<<std::endl; //", subrun: "<<storage->subrun_id()<<", "<<storage->last_subrun_id() <<std::endl ;


    _event++ ;

    // Now calculate the total POT + total numu neutrinos 
    auto ev_pot = storage->get_subrundata<potsummary>("generator"); 

    if( storage->subrun_id() != storage->last_subrun_id() )
      _tot_pot += ev_pot->totgoodpot ;

    auto ev_mctruth = storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) return false;

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();


    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    if( xyz[0] < _xmin || xyz[0] > _xmax || xyz[1] < _ymin || xyz[1] > _ymax || xyz[2] < _zmin || xyz[2] > _zmax ){
      return false;
      }

    _tot_event_in_AV ++; 

    auto parts = ev_mctruth->at(0).GetParticles();
    bool pi0 = false;
    bool mu  = false ;
    int n_pi0 = 0;
    int n_mes = 0;
    int n_lep = 0;
    int n_mu = 0;
      
    for ( auto const & p : parts ){
      
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){
        pi0 = true;
	n_pi0 += 1;
	}

      if( p.StatusCode() == 1 && p.PdgCode() == 13 ){
        n_mu += 1;
        mu = true;  
	}

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 211 || abs(p.PdgCode()) == 321 ||  p.PdgCode() == 130 
                                   || p.PdgCode() == 310 || abs(p.PdgCode()) == 311 ) )
        n_mes += 1;

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 11 || p.PdgCode() == -13))
        n_lep += 1;
      }

      if( n_mu == 1 && n_pi0 == 1 ){ //&& n_lep == 0 && n_mes == 0){
        //std::cout<<"**********************FOUND CCPI0!! "<<std::endl;
        _signal++;
	_event_list.emplace_back(_event-1);
        }

    return true;
  }

  bool CrossSection::finalize() {

     float FV = 216.35 * 193 * 1016.8 ;  // Fiducial volume 
     //float AV = 256.35 * 232. * 1036.8 ; // Active volume 

     float rho = 1.4 ; // g / cm3
     float avogadro = 6.022*pow(10,23);
     float g_per_mole = 39.948 ;
     float n_nucleon = 40 ;
     float n_targ = rho * FV / g_per_mole * avogadro * n_nucleon ;

     float flux = 8.72*pow(10,9) ; //From technote float(_n_numu) / FA / _tot_pot ;
     //std::cout<<"CCpi0 are "<<float(_signal)/(_tot_event_in_AV)*100<<"\% of BNB ("<<_signal<<"/"<<_tot_event_in_AV<<")"<<std::endl ;
     std::cout<<"CCpi0 are "<<float(_signal)/_event*100<<"\% of BNB ("<<_signal<<"/"<<_event<<")"<<std::endl ;

     _mean_e /= _n_numu ;
      
     //std::cout<<"Signal in FV: "<<_signal<<std::endl ;
     //std::cout<<"N targets: "<<n_targ <<std::endl ;
     //std::cout<<"Total POT: "<<_tot_pot<<", flux: "<<flux<<" numu/cm3/POT"<<std::endl;

     std::cout<<"Total POT: "<<_tot_pot<<std::endl;
     //std::cout<<"MC Cross section : "<<float(_signal) / n_targ / flux ; 

    std::cout<<_event_list.size()<<" in Event list :" <<std::endl ;
    for( auto const & e : _event_list ) std::cout<<e<<", ";


  
    return true;
  }

}
#endif
