#ifndef LARLITE_GETEVENTWEIGHTSFROMFULLSAMPLE_CXX
#define LARLITE_GETEVENTWEIGHTSFROMFULLSAMPLE_CXX

#include "GetEventWeightsFromFullSample.h"
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

  GetEventWeightsFromFullSample::GetEventWeightsFromFullSample() {
   
    _name                    = "GetEventWeightsFromFullSample";
    _fout                    = 0;
    _eventweight_label       = "" ;
    _n_sig = 0;

    _map_v.resize(10,std::multimap<float,float>());

    _events = 0;
    _event_no_dup = 0;


  }
  
  bool GetEventWeightsFromFullSample::initialize(){

      // Flux stuff
      _flux_label_v = {"expskin_FluxUnisim", "horncurrent_FluxUnisim", "kminus_PrimaryHadronNormalization", "kplus_PrimaryHadronFeynmanScaling", "kzero_PrimaryHadronSanfordWang", "nucleoninexsec_FluxUnisim", "nucleonqexsec_FluxUnisim", "nucleontotxsec_FluxUnisim", "piminus_PrimaryHadronSWCentralSplineVariation", "pioninexsec_FluxUnisim",  "pionqexsec_FluxUnisim",  "piontotxsec_FluxUnisim",  "piplus_PrimaryHadronSWCentralSplineVariation"};

      _unisim_label_v = {"expskin_FluxUnisim", "horncurrent_FluxUnisim", "nucleoninexsec_FluxUnisim", "nucleonqexsec_FluxUnisim", "nucleontotxsec_FluxUnisim", "pioninexsec_FluxUnisim",  "pionqexsec_FluxUnisim",  "piontotxsec_FluxUnisim"};

      _label_map = { {"kplus_PrimaryHadronFeynmanScaling", 1} , {"kminus_PrimaryHadronNormalization",2}, {"kzero_PrimaryHadronSanfordWang",3}, {"piplus_PrimaryHadronSWCentralSplineVariation",4}, {"piminus_PrimaryHadronSWCentralSplineVariation",5} } ;

      _funcs = (_flux_label_v.size() - _unisim_label_v.size() + 1);
      _t_weights_by_universe.resize(_flux_label_v.size() - _unisim_label_v.size() + 1);

      for ( int i = 0; i < _funcs ; i++ )
        _t_weights_by_universe.at(i).resize(1000,0) ;


      // GENIE stuff   
      _genie_label_v = {"AGKYpT","AGKYxF","DISAth","DISBth","DISCv1u","DISCv2u","FermiGasModelKf", "FermiGasModelSf","FormZone", "IntraNukeNabs", "IntraNukeNcex", "IntraNukeNel", "IntraNukeNinel", "IntraNukeNmfp", "IntraNukeNpi", "IntraNukePIabs", "IntraNukePIcex", "IntraNukePIel", "IntraNukePIinel", "IntraNukePImfp", "IntraNukePIpi", "NC", "NonResRvbarp1pi", "NonResRvbarp2pi", "NonResRvp1pi", "NonResRvp2pi", "ResDecayEta", "ResDecayGamma", "ResDecayTheta", "ccresAxial", "ccresVector", "cohMA", "cohR0", "ncelAxial", "ncelEta", "ncresAxial", "ncresVector", "qema", "qevec"};

      int genie_funcs = _genie_label_v.size() ;
     _all_evts_nominal = 0;
     _all_evts_m1.resize(genie_funcs,0) ; 
     _all_evts_p1.resize(genie_funcs,0) ; 

      //_t_weights_by_universe.resize(_genie_label_v.size());
      //_label_map={ {"genie_AGKYpT_Genie",0},{"genie_AGKYxF_Genie",1},{"genie_DISAth_Genie",2},{"genie_DISBth_Genie",3},{"genie_DISCv1u_Genie",4},
      //            {"genie_DISCv2u_Genie",5},{"genie_FermiGasModelKf_Genie",6},{"genie_FermiGasModelSf_Genie",7},{"genie_FormZone_Genie",8},
      //     {"genie_IntraNukeNabs_Genie",9},{"genie_IntraNukeNcex_Genie",10},{"genie_IntraNukeNel_Genie",11},{"genie_IntraNukeNinel_Genie",12},
	  // {"genie_IntraNukeNmfp_Genie",13},{"genie_IntraNukeNpi_Genie",14},{"genie_IntraNukePIabs_Genie",15},
	  // {"genie_IntraNukePIcex_Genie",16},{"genie_IntraNukePIel_Genie",17},{"genie_IntraNukePIinel_Genie",18},{"genie_IntraNukePImfp_Genie",19},
	  // {"genie_IntraNukePIpi_Genie",20},{"genie_NC_Genie",21},{"genie_NonResRvbarp1pi_Genie",22},
	  // {"genie_NonResRvbarp2pi_Genie",23},{"genie_NonResRvp1pi_Genie",24},{"genie_NonResRvp2pi_Genie",25},{"genie_ResDecayEta_Genie",26},
	  // {"genie_ResDecayGamma_Genie",27},{"genie_ResDecayTheta_Genie",28},{"genie_ccresAxial_Genie",29},
	  // {"genie_ccresVector_Genie",30},{"genie_cohMA_Genie",31},{"genie_cohR0_Genie",32},{"genie_ncelAxial_Genie",33},{"genie_ncelEta_Genie",34},
	  // {"genie_ncresAxial_Genie",35},{"genie_ncresVector_Genie",36},{"genie_qema_Genie",37},{"genie_qevec_Genie",38}};


    std::cout<<"EVENT WEIGHT LABELS, flux + genie : "<<_flux_label_v.size()<<", "<<genie_funcs<<std::endl ;

    return true;
   }

 
  bool GetEventWeightsFromFullSample::analyze(storage_manager* storage) {

      _events ++ ;

      auto r = storage->run_id() ;
      auto it = _map_v.at(r).find(storage->subrun_id());
      bool foundit = false;

      if( it != _map_v.at(r).end() ){
       while ( it->first == storage->subrun_id() ){  
         auto temp_event = it->second ; 
         if( temp_event == storage->event_id() )
           foundit = true;

         it++; 
         }   
        if ( !foundit)
         _map_v.at(r).emplace(storage->subrun_id(), storage->event_id() );
       
        else{
          //std::cout<<"Event: "<<_events<<std::endl;
          //std::cout<<"Duplicates "<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl;
          return false ;
        }   
       }   
      else 
         _map_v.at(r).emplace(storage->subrun_id(), storage->event_id() );

    auto ev_mctruth = storage->get_data<event_mctruth>("generator"); 
    auto ev_gen_wgt= storage->get_data<event_mceventweight>("genieeventweight"); 
    auto ev_flux_wgt= storage->get_data<event_mceventweight>("fluxeventweight"); 

    if(!ev_mctruth || !ev_mctruth->size() ){ 
       std::cout<<"No Truth..."<<std::endl ;
       return false;
     }

    if( !ev_gen_wgt || !ev_gen_wgt->size() ){
      std::cout<<"No GENIE event weights..." <<std::endl;
      return false;
     }

    if( !ev_flux_wgt || !ev_flux_wgt->size() ){
      std::cout<<"No FLUX event weights..." <<std::endl;
      return false;
     }

    _event_no_dup++ ;

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


    std::vector<float> unisim_tot_weight(1000,1);
    
    if ( _signal ) { // && _eventweight_label == "fluxeventweight" ) {
      auto flux_wgt  = ev_flux_wgt->at(0).GetWeights();

      for ( auto const & m : flux_wgt ) {
        std::cout<<"Size fo weights "<<m.first<<", "<<m.second.size()<<std::endl ;
	    if (m.first == "bnbcorrection_FluxHist" )
	      continue;
  	    // There should be 1000 weights
        for ( int jj = 0; jj < m.second.size(); jj++){
          
          // These are the 5 hadron values -- they may be correlated. Vary these 5 individually.
          if ( find(_unisim_label_v.begin(),_unisim_label_v.end(),m.first) != _unisim_label_v.end() ){
               unisim_tot_weight[jj] *= m.second.at(jj) ;
               //if ( jj < 10) std::cout<<"IS fluxunisim: "<<m.first<<", iteration: "<<jj<<", "<<m.second.at(jj)<<", "<<unisim_tot_weight[jj]<<std::endl;
          }
          else{
	      //std::cout<<"_label stuff: "<<m.first<<", "<<std::endl;
              _t_weights_by_universe[_label_map[m.first]][jj] += m.second.at(jj);
              //if ( jj < 10) std::cout<<"NOT fluxunisim: "<<m.first<<", iteration: "<<jj<<", "<<m.second.at(jj)<<", "<<_t_weights_by_universe[_label_map[m.first]][jj]<<std::endl;
	      }
        }
      }
      for ( int jj = 0; jj < 1000; jj++)
        _t_weights_by_universe[0][jj] += unisim_tot_weight[jj] ;

    }

    if ( _signal ){ //&& _eventweight_label == "genieeventweight" ) {

      _all_evts_nominal ++ ;

      int it = 0;
      auto gen_wgt  = ev_gen_wgt->at(0).GetWeights();

      // According to genie weight package, plus sig is stored first 
      for ( auto const & m : gen_wgt ) { 
         auto w_v = m.second ;
         //std::cout<<"Parameter: "<<m.first<<", "<<m.second.at(0)<<", "<<m.second.at(1)<<std::endl;
         _all_evts_p1[it] += (w_v.at(0)) ; 
         _all_evts_m1[it] += (w_v.at(1)) ;
         it++;
       }   

    }

    return true;
  }

  bool GetEventWeightsFromFullSample::finalize() {

    std::cout<<"Events: "<<_events <<", "<<_event_no_dup<<", "<<_n_sig<<std::endl ;


    if ( _print_output ) { // && _eventweight_label == "fluxeventweight" ){

      std::cout<<"{";

      for( int i = 0; i < _funcs; i++){
        std::cout<<"{";
        for( int j= 0; j < 1000; j++){
          if ( j != 999 )
            std::cout<<  _t_weights_by_universe[i][j]<<", ";
          else 
            std::cout<<  _t_weights_by_universe[i][j];
        }
        if ( i != _funcs - 1 ) std::cout<<"},";
        else std::cout<<"}";
      }
      std::cout<<"}";

      std::cout<<std::endl ;
    }

    if ( _print_output ) { //&& _eventweight_label == "genieeventweight" ){

      std::cout<<"All GENIE events: "<<_all_evts_nominal<<std::endl ;
      for( int i = 0 ; i < _all_evts_m1.size(); i++) {
        std::cout<<"\nFunction: "<<_genie_label_v[i]<<std::endl ;
        std::cout<<"All events (-3sig): "<<_all_evts_m1.at(i)<<std::endl ;
        std::cout<<"All events (+3sig): "<<_all_evts_p1.at(i)<<std::endl ;
      }

      std::cout<<"All events: "<<_all_evts_nominal<<std::endl ;

      for( int i = 0 ; i < _all_evts_m1.size(); i++)
        std::cout<<_all_evts_m1[i]<<", " ;

      std::cout<<std::endl ;
      for( int i = 0 ; i < _all_evts_p1.size(); i++)
        std::cout<<_all_evts_p1[i]<<", " ;

      std::cout<<std::endl ;

    }

    return true;
  }

}
#endif
