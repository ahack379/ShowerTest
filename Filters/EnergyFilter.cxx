#ifndef LARLITE_ENERGYFILTER_CXX
#define LARLITE_ENERGYFILTER_CXX

#include "EnergyFilter.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool EnergyFilter::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool EnergyFilter::analyze(storage_manager* storage) {
    
    auto ev_truth = storage->get_data<event_mctruth>("generator");
    auto nu_truth = ev_truth->at(0).GetNeutrino();
    auto nu_part_traj = nu_truth.Nu().Trajectory();

    auto vtxX = nu_part_traj.at( nu_part_traj.size() - 1).X(); 
    auto vtxY = nu_part_traj.at( nu_part_traj.size() - 1).Y(); 
    auto vtxZ = nu_part_traj.at( nu_part_traj.size() - 1).Z(); 

   // auto ev_vtx = storage->get_data<event_vertex>("pandoraNu");

    if( vtxX < 20 || vtxX > 236.5 || vtxY < -96.5 || vtxY > 96.5 || vtxZ < 10 || vtxZ > 1026.8 )
      return false ;

    auto parts = ev_truth->at(0).GetParticles();
    bool pi0 = false;

    float energy = 0;

    for ( auto const & p : parts ){

      //std::cout<<"Each pdg: "<<p.PdgCode()<<", "<<p.Process()<<std::endl ;
    
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){  
        pi0 = true;
        energy = p.Trajectory().at(0).E();
        //std::cout<<"Energy: "<<energy<<std::endl ;
        if( energy > .200 ){
          std::cout<<"REMOVE: "<<energy<<std::endl;
          return false ;
          }   
        }
      }
  
    return true;
  }

  bool EnergyFilter::finalize() {

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
