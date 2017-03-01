#ifndef LARLITE_SIGNALEFF_CXX
#define LARLITE_SIGNALEFF_CXX

#include "SignalEff.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool SignalEff::initialize() {

    _event = 0 ;
    _sel_good = 0 ;
    _sel_misid = 0 ;
    _sel_tot = 0 ;

    return true;
  }
  
  bool SignalEff::analyze(storage_manager* storage) {

    // This module runs on Sel2 output 
    // How many of the pi0 events in this set do we get?
    // PURPOSE OF MODULE: Access purity and efficiency of selectionII using MCC8 adjustments
    auto ev_truth = storage->get_data<event_mctruth>("generator");
    auto & truth = ev_truth->at(0);
    auto & nu  = truth.GetNeutrino();
    
    auto t = nu.InteractionType();

    bool isSel = IsSelected(_sel_ev_list, _event);
     
    if( isSel && (t == 1004 || t == 1011 || t == 1080 || t == 1086 || t == 1090) ){ 
      _sel_good++ ;
      std::cout<<"Ineraction Type: "<<nu.InteractionType()<<std::endl ;
       }

    if( isSel && (t != 1004 && t != 1011 && t != 1080 && t != 1086 && t != 1090) )
      _sel_misid++;

    if( t == 1004 || t == 1011 || t == 1080 || t == 1086 || t == 1090 )
      _sel_tot++;
    
    _event++;

    return true;
  }

  bool SignalEff::finalize() {

    std::cout<<"Efficiency: "<<float(_sel_good)/_sel_tot*100<<"% ("<<_sel_good<<"/"<<_sel_tot<<")"<<std::endl ;
    std::cout<<"MisID Rate: "<<float(_sel_misid)/_sel_tot*100<<"% ("<<_sel_misid<<"/"<<_sel_tot<<")"<<std::endl;

    /// Need to normalize this number
    std::cout<<"Cosmic Background: "<<_n_cosbkgd_evts<<std::endl ;
  
    return true;
  }

  bool SignalEff::IsSelected(std::vector<int> sel_list, const int & ev){

    for( size_t i = 0; i < sel_list.size(); i++ ){
      if( ev == sel_list[i] )
        return true;
      }
     return false; 
  }

}
#endif
