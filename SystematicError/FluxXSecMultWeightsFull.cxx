#ifndef LARLITE_FLUXXSECMULTWEIGHTSFULL_CXX
#define LARLITE_FLUXXSECMULTWEIGHTSFULL_CXX

#include "FluxXSecMultWeightsFull.h"
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

#include "TH1D.h"
#include "TDirectory.h"



namespace larlite {

  FluxXSecMultWeightsFull::FluxXSecMultWeightsFull() {
   
    _name                    = "FluxXSecMultWeightsFull";
    _fout                    = 0;

    _n_sig = 0;

   _genie_label_v = {"expskin_FluxUnisim", "horncurrent_FluxUnisim", "kminus_PrimaryHadronNormalization", "kplus_PrimaryHadronFeynmanScaling", "kzero_PrimaryHadronSanfordWang", "nucleoninexsec_FluxUnisim", "nucleonqexsec_FluxUnisim", "nucleontotxsec_FluxUnisim", "piminus_PrimaryHadronSWCentralSplineVariation", "pioninexsec_FluxUnisim",  "pionqexsec_FluxUnisim",  "piontotxsec_FluxUnisim",  "piplus_PrimaryHadronSWCentralSplineVariation"};

   _unisim_label_v = {"expskin_FluxUnisim", "horncurrent_FluxUnisim", "nucleoninexsec_FluxUnisim", "nucleonqexsec_FluxUnisim", "nucleontotxsec_FluxUnisim", "pioninexsec_FluxUnisim",  "pionqexsec_FluxUnisim",  "piontotxsec_FluxUnisim"};

    _label_map = { {"kplus_PrimaryHadronFeynmanScaling", 1} , {"kminus_PrimaryHadronNormalization",2}, {"kzero_PrimaryHadronSanfordWang",3}, {"piplus_PrimaryHadronSWCentralSplineVariation",4}, {"piminus_PrimaryHadronSWCentralSplineVariation",5} } ;


   //_unisim_label_v = {"SkinEffect","HornCurrent","NucleonInXsec","NucleonQEXsec","NucleonTotXsec","piInelasticXsec","piInelasticXsec","piQEXsec","piTotalXsec"};

  int funcs = (_genie_label_v.size() - _unisim_label_v.size() + 1);

  _t_weights_by_universe.resize(_genie_label_v.size() - _unisim_label_v.size() + 1);

  for ( int i = 0; i < funcs ; i++ )
      _t_weights_by_universe.at(i).resize(1000,0) ;
  

    std::cout<<"GENIE LABLES: "<<_genie_label_v.size()<<std::endl ;
    
    fGeometry = nullptr;
    _events = 0;

  }
  
  bool FluxXSecMultWeightsFull::initialize(){

    return true;
   }

 
  bool FluxXSecMultWeightsFull::analyze(storage_manager* storage) {

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

   int n_pi0 = 0;
   int n_mu = 0;

   auto parts = ev_mctruth->at(0).GetParticles();
   for ( auto const & p : parts ){
   
     if( p.StatusCode() == 1 && p.PdgCode() == 111 )
       n_pi0 += 1;

     if( p.StatusCode() == 1 && p.PdgCode() == 13 )
       n_mu += 1;
    }

    if ( n_pi0 == 1 && n_mu == 1 && infv ){
      _signal = true ;
      _n_sig ++ ;
    }
    else
      _signal = false;

    auto wgt  = ev_wgt->at(0).GetWeights();

    std::vector<float> unisim_tot_weight(1000,1);
    
    if ( _signal ) {
      for ( auto const & m : wgt ) {
        //std::cout<<"Size fo weights "<<m.first<<", "<<m.second.size()<<std::endl ;
        for ( int jj = 0; jj < m.second.size(); jj++){
          // THese are the 5 hadron values -- they may be correlated. Vary these 5 individually.
          if ( find(_unisim_label_v.begin(),_unisim_label_v.end(),m.first) != _unisim_label_v.end() )
               unisim_tot_weight[jj] *= m.second.at(jj) ;
          else
              _t_weights_by_universe[_label_map[m.first]][jj] += m.second.at(jj);
        }
      }
      for ( int jj = 0; jj < 1000; jj++)
        _t_weights_by_universe[0][jj] += unisim_tot_weight[jj] ;

    }
      

    return true;
  }

  bool FluxXSecMultWeightsFull::finalize() {

    std::cout<<"Events: "<<_events <<", "<<_n_sig<<std::endl ;

    int funcs = _genie_label_v.size() - _unisim_label_v.size() + 1 ; 

    std::cout<<"{";

    for( int i = 0; i < funcs; i++){
      std::cout<<"{";
      for( int j= 0; j < 1000; j++){
        if ( j != 999 )
          std::cout<<  _t_weights_by_universe[i][j]<<", ";
        else 
          std::cout<<  _t_weights_by_universe[i][j];
      }
      if ( i != funcs - 1 ) std::cout<<"},";
      else std::cout<<"}";
    }
    std::cout<<"}";

    std::cout<<std::endl ;

    return true;
  }

}
#endif
