#ifndef LARLITE_SPLITDATAOddEven_CXX
#define LARLITE_SPLITDATAOddEven_CXX

#include "SplitDataOddEven.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/potsummary.h"

namespace larlite {

  bool SplitDataOddEven::initialize() {

    if(!_tree){
      _tree = new TTree("tree","tree");
      _tree->Branch("x_even",&_x_even,"x_even/I");
      _tree->Branch("x_odd",&_x_odd,"x_odd/I");
      _tree->Branch("x_even_POT",&_x_even_POT,"x_even_POT/F");
      _tree->Branch("x_odd_POT",&_x_odd_POT,"x_odd_POT/F");
    }

    _x_even = 0;
    _x_odd = 0;

    _x_even_POT = 0;
    _x_odd_POT = 0;

    return true;
  }
  
  bool SplitDataOddEven::analyze(storage_manager* storage) {
  
    auto event_id = storage->event_id() ;

    bool IsEven = event_id % 2 == 0 ? 1 : 0 ;

    //std::cout<<"event id : "<<event_id <<", "<<IsEven <<std::endl;

    // Now calculate the total POT + total numu neutrinos 
    auto ev_pot = storage->get_subrundata<potsummary>("generator"); 

    if ( IsEven ){
      _x_even ++ ;
      if( storage->subrun_id() != storage->last_subrun_id() )
        _x_even_POT += ev_pot->totgoodpot ;
    }
    else{
      _x_odd ++ ;
      if( storage->subrun_id() != storage->last_subrun_id() )
        _x_odd_POT += ev_pot->totgoodpot ;
    }

    return true;
  }

  bool SplitDataOddEven::finalize() {
    
    std::cout<<"Odd, Even: "<<_x_odd<<", "<<_x_even<<std::endl;
    std::cout<<"POT Odd, Even: "<<_x_odd_POT<<", "<<_x_even_POT<<std::endl;

    _tree->Fill();

    if(_fout) { _fout->cd(); _tree->Write(); }
  
    return true;
  }

}
#endif
