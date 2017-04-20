#ifndef LARLITE_XSECERROR_CXX
#define LARLITE_XSECERROR_CXX

#include "XSecError.h"

namespace larlite {

  bool XSecError::initialize() {

    if(!_tree){
    _tree = new TTree("tree","");
    _tree->Branch("test",&_test,"test/F");
    }

    _hist = new TH1F("hist","hist",20,0,50);
    _hist2 = new TH1F("hist2","hist2",20,0,50);
    _hist3 = new TH1F("hist3","hist3",20,0,50);
    _histw = new TH1F("histw","histw",20,0,50);

    _test = 0;
    return true;
  }
  
  bool XSecError::analyze(storage_manager* storage) {
  
    _test++;

    for ( int i = 0; i < 10; i++){

      _hist->Fill(_test);
      _hist3->Fill(_test,4);
        if(i > 6 ){
        _hist2->Fill(_test);
        _histw->Fill(_test,4);
        }
      }
  
    return true;
  }

  bool XSecError::finalize() {

    if(_fout) { _fout->cd(); _hist->Write(); _histw->Write(); _hist2->Write(); _hist3->Write(); }

    std::cout<<"Unweighted INTEGRALS: "<<_hist->Integral()<<", "<<_hist2->Integral()<<std::endl ;
    std::cout<<"Weighted INTEGRALS: "<<_hist3->Integral()<<", "<<_histw->Integral()<<std::endl ;
    
    std::cout<<"Unweighted RATIOS : "<<float(_hist2->Integral())/_hist->Integral()<<std::endl ;
    std::cout<<"Weighted RATIOS :"<<float(_histw->Integral())/_hist3->Integral()<<std::endl ;
  
    return true;
  }

}
#endif
