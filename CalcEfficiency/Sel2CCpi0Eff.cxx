#ifndef LARLITE_SEL2CCPI0EFF_CXX
#define LARLITE_SEL2CCPI0EFF_CXX

#include "Sel2CCpi0Eff.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

namespace larlite {

  bool Sel2CCpi0Eff::initialize() {

  std::cout<<"Entering initialize "<<std::endl ;

    _events = -1 ;
    _signal = 0;

    if ( ! _cut_tree ){
      _cut_tree = new TTree("cuttree","cuttree");
      _cut_tree->Branch("_energy",&_energy,"_energy/F");
    }


    _event_list.clear();

    return true;
  }
  
  bool Sel2CCpi0Eff::analyze(storage_manager* storage) {

    _events++;

    std::cout<<"Event is : "<<_events <<std::endl ;

    // This is a calculation of CCpi0 in the BNB file 
    // Can also be run on sel2 filter to calculate the
    // eff of sel2 filter for OUR signal
    auto ev_truth = storage->get_data<event_mctruth>("generator");
    if(!ev_truth || !ev_truth->size() ) {std::cout<<"*********************"<<std::endl; return false; } 

    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    auto ev_mct = storage->get_data<event_mctrack>("mcreco");

    auto & truth = ev_truth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();


    bool pi0 = false ;

    for ( auto const & s : *ev_mcs ){

       //std::cout<<"Shower PDG: "<<s.PdgCode()<<std::endl ;
       if ( s.MotherPdgCode() == 111 ){

         auto st = s.Start() ;
         auto dist = sqrt( pow(xyz[0] - st.X(),2) + pow(xyz[1] - st.Y(),2) +pow(xyz[2] - st.Z(),2) ); 

          if( dist < 0.4 ) {
	    if( _CCNC )
	      _signal++;
	    std::cout<<"Found a pi0 with "<<!nu.CCNC()<<" leptons!!!! "<<std::endl ;
	    pi0 = true;
	    break;
           }
         }
       }

    if ( !pi0 )
      return false ;

    if( !_CCNC ) {

    for ( auto const & t : *ev_mct ){

       //std::cout<<"Track PDG: "<<t.PdgCode()<<std::endl ;
       if ( t.PdgCode() == 13 || t.PdgCode() == -13 ){

         auto st = t.Start() ;
         auto dist = sqrt( pow(xyz[0] - st.X(),2) + pow(xyz[1] - st.Y(),2) +pow(xyz[2] - st.Z(),2) ); 

          if( dist < 0.4 ) {

            _signal++ ;
            _event_list.emplace_back(_events);
            std::cout<<"Foudn a CCpi0!!!!!!!!!!!"<<std::endl ;

          bool found = false;
	  float min_e = 1e9;

          // This chunk of code finds the min energy of a pi0 shower and stores it in a tree for cuts based analysis

        //  for ( auto const & s : *ev_mcs ){
        //    if ( s.MotherPdgCode() == 111 ){
        //      auto st = s.Start() ;
        //      auto dist = sqrt( pow(xyz[0] - st.X(),2) + pow(xyz[1] - st.Y(),2) +pow(xyz[2] - st.Z(),2) ); 
        //        if( dist < 0.4 ) {
	//	  std::cout<<"Energy is : "<<s.DetProfile().E()<<std::endl ;

	//	  // If we enter this if block, we've found the second gamma 
        //          if( found == true ){

        //            _energy = s.DetProfile().E();

	//	    if( min_e < _energy )
	//	      _energy = min_e;

        //            if( _energy <= 30 ){
        //              _signal++ ;
        //              _event_list.emplace_back(_events);
	//	      }

        //            _cut_tree->Fill();
        //            break;
        //            }
	//	  min_e = s.DetProfile().E();
	//	   
        //          found = true;
        //          }
        //        }
        //      }
                  return true;
            }


          }
        }
      }


  }

  bool Sel2CCpi0Eff::finalize() {

    std::cout<<"CCpi0 are "<<float(_signal)/(_events+1)*100<<"\% of BNB ("<<_signal<<"/"<<_events+1<<")"<<std::endl ;
    if(_fout){
      _fout->cd();
      _cut_tree->Write();
      }

    std::cout<<_event_list.size()<<" in Event list :" <<std::endl ;
    for( auto const & e : _event_list )
      std::cout<<e<<", ";
  
    return true;
  }

}
#endif
