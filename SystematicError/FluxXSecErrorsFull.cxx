#ifndef LARLITE_FLUXXSECERRORSFULL_CXX
#define LARLITE_FLUXXSECERRORSFULL_CXX

#include "FluxXSecErrorsFull.h"
#include "DataFormat/track.h"
#include "DataFormat/opflash.h"
#include "DataFormat/vertex.h"
#include "DataFormat/calorimetry.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/hit.h"
#include "DataFormat/shower.h"
#include "DataFormat/potsummary.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mceventweight.h"
#include "DataFormat/gtruth.h"

// Here we'll do the flux renormalization 
// For each energy, we weight wiht fXsec 


namespace larlite {

  FluxXSecErrorsFull::FluxXSecErrorsFull() {
   
    _name                    = "FluxXSecErrorsFull";
    _fout                    = 0;

   _genie_label_v = {"expskin_FluxUnisim", "horncurrent_FluxUnisim", "kminus_PrimaryHadronNormalization", "kplus_PrimaryHadronFeynmanScaling", "kzero_PrimaryHadronSanfordWang", "nucleoninexsec_FluxUnisim", "nucleonqexsec_FluxUnisim", "nucleontotxsec_FluxUnisim", "piminus_PrimaryHadronSWCentralSplineVariation", "pioninexsec_FluxUnisim",  "pionqexsec_FluxUnisim",  "piontotxsec_FluxUnisim",  "piplus_PrimaryHadronSWCentralSplineVariation"};

    fGeometry = nullptr;
    _tree = nullptr;

  }
  
  bool FluxXSecErrorsFull::initialize(){

    fGeometry = larutil::Geometry::GetME();
    _tot_pot = 0. ;

    int funcs = _genie_label_v.size() ; 

    _all_evts_nominal = 0;
    _all_evts_m1.resize(funcs,0) ; 
    _all_evts_p1.resize(funcs,0) ; 

    if( !_tree) {
       _tree = new TTree("flux_tree","");
       _tree->Branch("cv",&_cv,"cv/F"); 
       _tree->Branch("up","std::vector<float>",&_up); 
       _tree->Branch("down","std::vector<float>",&_down); 
       _tree->Branch("gtruth_wgt",&_gtruth_wgt,"gtruth_wgt/F"); 
       _tree->Branch("nu_pdg",&_nu_pdg,"nu_pdg/I");
       _tree->Branch("signal",&_signal,"signal/B");
     }

    _cv = -1;
    _up.resize(funcs,-1);
    _down.resize(funcs,-1);

    _t_weights_by_universe.resize(funcs);

    for ( int i = 0; i < funcs ; i++ )
      _t_weights_by_universe.at(i).resize(1000,0) ;

    _file.open("list.txt",std::ios_base::in);
 
    return true;
   }

 
  bool FluxXSecErrorsFull::analyze(storage_manager* storage) {

    auto ev_mctruth = storage->get_data<event_mctruth>("generator"); 
    auto ev_wgt= storage->get_data<event_mceventweight>("fluxeventweight"); 
    auto ev_gtruth = storage->get_data<event_gtruth>("generator"); 

    _gtruth_wgt = ev_gtruth->at(0).fXsec ;

     if(!ev_mctruth || !ev_mctruth->size() ){ 
       std::cout<<"No Truth..."<<std::endl ;
       return false;
     }

    if( !ev_wgt || !ev_wgt->size() ){
      std::cout<<"No event weights..." <<std::endl;
      return false;
     }

    auto wgt  = ev_wgt->at(0).GetWeights();

    auto nu  = ev_mctruth->at(0).GetNeutrino();
    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();
    auto nue = traj.at(traj.size() - 1).E();
    _cv = nue ;
    _nu_pdg = nu.Nu().PdgCode();

    bool infv = true;
    if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
       infv = false ;

    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0;
    int n_mu = 0;
    
    for ( auto const & p : parts ){
    
      if( p.StatusCode() == 1 && p.PdgCode() == 111 )
        n_pi0 += 1;

      if( p.StatusCode() == 1 && p.PdgCode() == 13 )
        n_mu += 1;
    }

    if ( n_pi0 == 1 && n_mu == 1 && infv )
      _signal = true ;
    else 
      _signal = false;


    std::vector<float> mean_v(13,0) ;
    int kk = 0;

    for ( auto const & m : wgt ) {
      for ( int jj = 0; jj < m.second.size(); jj++){

        if (_signal)
          _t_weights_by_universe[kk][jj] += m.second.at(jj);

        mean_v[kk] += (m.second.at(jj)/1000);
      }
      kk++;
    }

    std::vector<float> sigma_v(13,0) ;
    kk = 0;

    for ( auto const & m : wgt ) {
      for ( int jj = 0; jj < m.second.size(); jj++){
        sigma_v[kk] += ( m.second.at(jj) - mean_v.at(kk) ) * ( m.second.at(jj) - mean_v.at(kk));
      }

      sigma_v[kk] /= 1000;

      _up[kk] = 1 + sqrt( sigma_v[kk] ) ;
      _down[kk] = 1 - sqrt( sigma_v[kk] ) ;

      kk ++;
    }
    _tree->Fill();


    return true;
  }

  bool FluxXSecErrorsFull::finalize() {

    auto funcs = _genie_label_v.size();

    std::cout<<"{";

    for( int i = 0; i < funcs; i++){
      std::cout<<"{";
      for( int j= 0; j < 1000; j++){
        if ( j != 999 )
          std::cout<<  _t_weights_by_universe[i][j]<<", ";
        else 
          std::cout<<  _t_weights_by_universe[i][j];
      }
      std::cout<<"}";
    }

    std::cout<<std::endl ;

    if ( _fout ){
      _fout->cd();
      _tree->Write();
    
    }

    return true;
  }

}
#endif
