#ifndef LARLITE_GETINTERACTIONINFO_CXX
#define LARLITE_GETINTERACTIONINFO_CXX

#include "GetInteractionInfo.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

namespace larlite {

  bool GetInteractionInfo::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool GetInteractionInfo::analyze(storage_manager* storage) {

    // eff of sel2 filter for OUR signal

    std::cout<<"\n\nNEW EVENT! "<<std::endl;
    auto ev_truth = storage->get_data<event_mctruth>("generator");
    if(!ev_truth || !ev_truth->size() ) return false;  
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    auto ev_mct = storage->get_data<event_mctrack>("mcreco");

    auto & truth = ev_truth->at(0);
    auto & nu  = truth.GetNeutrino();

    auto & traj0 = nu.Nu().Trajectory();

    auto const & E  = traj0.at(traj0.size() - 1).E();
    std::cout<<"Neutrino of energy "<<E<<" is of interact type "<<nu.InteractionType()<<std::endl ;

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    for ( auto const & s : *ev_mcs ){
       std::cout<<"Shower Pdgcode : "<<s.PdgCode()<<", "<<s.DetProfile().E()<<", mother pdg: "<<s.MotherPdgCode()<<std::endl;

       if ( s.MotherPdgCode() == 111 ){
         auto st = s.Start() ;
         auto dist = sqrt( pow(xyz[0] - st.X(),2) + pow(xyz[1] - st.Y(),2) +pow(xyz[2] - st.Z(),2) );
	 std::cout<<"Found a pi0...."<<dist<<std::endl;
         if( dist < 0.4 ) {
            std::cout<<"...from the origin with "<<!nu.CCNC()<<" leptons!!!! "<<std::endl ;
            //break;
           }
         }
       }

    for ( auto const & t : *ev_mct ){
       if ( t.size() > 5 )
         std::cout<<"Track Pdgcode : "<<t.PdgCode()<<std::endl;
	 }


    auto parts = ev_truth->at(0).GetParticles();
    for ( auto const & p : parts ){
      
      if( p.StatusCode() == 1){
         std::cout<<"Particle PDG "<<p.PdgCode()<<std::endl ; //", and Parent: "<<p.Mother()<<std::endl ;
	 }

      }

    return true;
  }

  bool GetInteractionInfo::finalize() {

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
