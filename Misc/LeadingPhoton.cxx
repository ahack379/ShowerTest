#ifndef LARLITE_LEADINGPHOTON_CXX
#define LARLITE_LEADINGPHOTON_CXX

#include "LeadingPhoton.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/shower.h"

namespace larlite {

  bool LeadingPhoton::initialize() {    

    _event = -1; 
    _one_shower_events = 0;

    if( !_tree ) {
      _tree = new TTree("tree","tree");
      _tree->Branch("dot_prod",&_dot_prod,"dot_prod/F");
      _tree->Branch("reco_e",&_reco_e,"reco_e/F");
      _tree->Branch("mc_e",&_mc_e,"mc_e/F");
    }

    return true;
  }
  
  bool LeadingPhoton::analyze(storage_manager* storage) {

    _event++ ;

    auto ev_s = storage->get_data<event_shower>("showerreco");
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");

    if( !ev_s || ev_s->size() != 1 ) {
      std::cout<<"Not enough reco'd showers..." <<std::endl; 
      return false;  
    }

    if( !ev_mcs || !ev_mcs->size() ) {
      std::cout<<"Not enough reco'd mcshowers..." <<std::endl; 
      return false;  
    }

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

    std::vector<int> shr_ids;
    float max_mce = -10;
    int max_it = -1;

    _one_shower_events ++ ;
   
    for ( int si = 0; si < ev_mcs->size(); si++){ 

      auto s = ev_mcs->at(si);

      if( s.PdgCode() != 22 ) 
        continue; 
    
      auto st = s.Start();
      auto dist = sqrt( pow(st.X() - xyz[0],2) + pow(st.Y() - xyz[1],2) + pow(st.Z() - xyz[2],2) );
    
      if ( dist < 0.0001 ){

        shr_ids.emplace_back(si) ;

        if ( s.DetProfile().E() > max_mce ){
          max_mce = s.Start().E() ;
          max_it = si ;
        }
      }
    }   
    
    if( shr_ids.size() ){

      auto mcs_i = ev_mcs->at(max_it) ;
      auto recos_i = ev_s->at(0) ;

      auto mag_mcs = sqrt( mcs_i.Start().Px()*mcs_i.Start().Px()+mcs_i.Start().Py()*mcs_i.Start().Py()+mcs_i.Start().Pz()*mcs_i.Start().Pz() );
      auto mag_reco = sqrt( pow(recos_i.Direction().Px(),2) + pow(recos_i.Direction().Py(),2) + pow(recos_i.Direction().Pz(),2) );

      _dot_prod = mcs_i.Start().Px() * recos_i.Direction().Px() +
                 mcs_i.Start().Py() * recos_i.Direction().Py() +
                 mcs_i.Start().Pz() * recos_i.Direction().Pz() ;

      _dot_prod /= ( mag_mcs * mag_reco );

      _reco_e = recos_i.Energy();
      _mc_e = mcs_i.Start().E() ;

      //std::cout<<"DOT is: "<<_dot_prod <<std::endl ;

      _tree->Fill();
    }


    return true;
  }

  bool LeadingPhoton::finalize() {

    //std::cout<<"\n\n"<<_event_list.size()<<" in Event list :" <<std::endl ;
    //for( auto const & e : _event_list) std::cout<<e<<", ";

    std::cout<<"One shower events: "<<_one_shower_events <<std::endl ;

    if ( _fout ){
      _fout->cd();
      _tree->Write();
    }
  
    return true;
  }

}
#endif
