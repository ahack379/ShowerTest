#ifndef LARLITE_CALCCHARGE_CXX
#define LARLITE_CALCCHARGE_CXX

#include "CalcCharge.h"
#include "DataFormat/hit.h"

namespace larlite {

  bool CalcCharge::initialize() {

   if( !_tree ){
     _tree = new TTree("tree","tree"); 
     _tree->Branch("sum_charge00",&_sum_charge00,"sum_charge00/F"); 
     _tree->Branch("sum_charge02",&_sum_charge02,"sum_charge02/F"); 
    }

    return true;
  }
  
  bool CalcCharge::analyze(storage_manager* storage) {
  
    auto ev_hit00 = storage->get_data<event_hit>("hit00");
    auto ev_hit02 = storage->get_data<event_hit>("hit02");

    if ( !ev_hit00 || !ev_hit00->size() ) return false; 

    if ( !ev_hit02 || !ev_hit02->size() ) return false; 

    _sum_charge00 = 0.;
    _sum_charge02 = 0.;

    for ( auto const & h : *ev_hit00 )
      _sum_charge00 += h.Integral(); 

    for ( auto const & h : *ev_hit02 )
      _sum_charge02 += h.Integral(); 

    _tree->Fill();
  
    return true;
  }

  bool CalcCharge::finalize() {

    if(_fout){
      _fout->cd(); 
      _tree->Write(); 
      }
  
    return true;
  }

}
#endif
