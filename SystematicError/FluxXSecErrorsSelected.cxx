#ifndef LARLITE_FLUXXSECERRORSSELECTED_CXX
#define LARLITE_FLUXXSECERRORSSELECTED_CXX

#include "FluxXSecErrorsSelected.h"
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

// HitRatio
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

namespace larlite {

  FluxXSecErrorsSelected::FluxXSecErrorsSelected() {
   
    _name                    = "FluxXSecErrorsSelected";
    _fout                    = 0;
    _tree = 0;

   _genie_label_v = {"expskin_FluxUnisim", "horncurrent_FluxUnisim", "kminus_PrimaryHadronNormalization", "kplus_PrimaryHadronFeynmanScaling", "kzero_PrimaryHadronSanfordWang", "nucleoninexsec_FluxUnisim", "nucleonqexsec_FluxUnisim", "nucleontotxsec_FluxUnisim", "piminus_PrimaryHadronSWCentralSplineVariation", "pioninexsec_FluxUnisim",  "pionqexsec_FluxUnisim",  "piontotxsec_FluxUnisim",  "piplus_PrimaryHadronSWCentralSplineVariation"};


     std::cout<<"GENIE LABLES: "<<_genie_label_v.size()<<std::endl ;
    
    fGeometry = nullptr;
    _event_producer = "genieeventweight" ;

    _events = 0;
  }
  
  bool FluxXSecErrorsSelected::initialize(){

    fGeometry = larutil::Geometry::GetME();
    _tot_pot = 0. ;
    int funcs = _genie_label_v.size() ; // 13 total for flux

    _sel_evts_nominal = 0;
    _sel_evts_m1.resize(funcs,0) ; 
    _sel_evts_p1.resize(funcs,0) ; 
    _sel_evts_nom.resize(funcs,0) ; 

    if( !_tree) {
       _tree = new TTree("flux_tree","");
       _tree->Branch("cv",&_cv,"cv/F"); 
       _tree->Branch("up","std::vector<float>",&_up); 
       _tree->Branch("down","std::vector<float>",&_down); 
     }

    _cv = -1;
    _up.resize(funcs,-1);
    _down.resize(funcs,-1);

    return true;
   }

 
  bool FluxXSecErrorsSelected::analyze(storage_manager* storage) {

    auto ev_mctruth = storage->get_data<event_mctruth>("generator"); 
    auto ev_wgt= storage->get_data<event_mceventweight>("fluxeventweight"); 

    if(!ev_mctruth || !ev_mctruth->size() ){ 
       std::cout<<"No Truth..."<<std::endl ;
       return false;
     }

    if( !ev_wgt || !ev_wgt->size() ){
      std::cout<<"No event weights..." <<std::endl;
      return false;
     }

    _events ++ ;

    auto nu  = ev_mctruth->at(0).GetNeutrino();
    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();
    auto nue = traj.at(0).E();
    _cv = nue ;

    bool infv = true;

   if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
      infv = false; 

    auto parts = ev_mctruth->at(0).GetParticles();

    std::vector<float> mean_v(13,0) ;
    int kk = 0;

    auto wgt  = ev_wgt->at(0).GetWeights();

    for ( auto const & m : wgt ) {
      //std::cout<<"Size fo weights: "<<m.second.size()<<std::endl ;
      for ( int jj = 0; jj < m.second.size(); jj++){
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
     //if ( kk == 12 ){
     //  std::cout<<"sig : "<<_up[kk]<<", "<<_down[kk]<<std::endl ;
     //  std::cout<<"Mean, sigma: "<< sigma_v[kk]<<", "<<mean_v[kk]<<std::endl ;
     //}

      kk ++;
    }
    _tree->Fill();

    _sel_evts_nominal ++ ;

      kk = 0 ;
      for ( auto const & m : wgt ) {
        
        _sel_evts_p1[kk] += ( (mean_v[kk] + sqrt(sigma_v[kk])) );
        _sel_evts_m1[kk] += ( (mean_v[kk] - sqrt(sigma_v[kk])) );
        _sel_evts_nom[kk] += ( mean_v[kk] );
        kk ++;
      }

    return true;
  }

  bool FluxXSecErrorsSelected::finalize() {

    std::cout<<"All events: "<<_sel_evts_nominal<<std::endl ;

    std::cout<<"Events: "<<_events <<std::endl ;

    for( int i = 0 ; i < _sel_evts_nom.size(); i++) 
      std::cout<<_sel_evts_nom[i]<<", " ;

    std::cout<<std::endl;

    for( int i = 0 ; i < _sel_evts_m1.size(); i++) 
      std::cout<<_sel_evts_m1[i]<<", " ;

    std::cout<<std::endl ;
    for( int i = 0 ; i < _sel_evts_p1.size(); i++) 
      std::cout<<_sel_evts_p1[i]<<", " ;

    std::cout<<std::endl ;

   if (_fout) {
    _fout->cd();
    _tree->Write();
   }

    return true;
  }

}
#endif
