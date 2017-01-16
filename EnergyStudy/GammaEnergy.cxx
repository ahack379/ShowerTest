#ifndef LARLITE_GAMMAENERGY_CXX
#define LARLITE_GAMMAENERGY_CXX

#include "GammaEnergy.h"
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/cluster.h"
#include "DataFormat/shower.h"
#include "DataFormat/event_ass.h"
#include "LArUtil/DetectorProperties.h"

namespace larlite {

  bool GammaEnergy::initialize() {

    if(!_tree){
      _tree = new TTree("tree","tree");
      _tree->Branch("hit_reco_e",&_hit_reco_e,"hit_reco_e/F"); 
      _tree->Branch("reco_e",&_reco_e,"reco_e/F"); 
      _tree->Branch("mc_e",&_mc_e,"mc_e/F"); 
      _tree->Branch("sum",&_sum,"sum/F"); 
      _tree->Branch("sum_adj",&_sum_adj,"sum_adj/F"); 
      _tree->Branch("sum_reco_int",&_sum_reco_int,"sum_reco_int/F"); 
     }

    _event = 0;

    return true;
  }
  
  bool GammaEnergy::analyze(storage_manager* storage) {
  
    std::cout<<"\nEvent is : "<<_event <<std::endl ;
    _event ++ ;

    auto ev_hit = storage->get_data<event_hit>("gaushit");
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    auto ev_s = storage->get_data<event_shower>("showerreco");
    auto ev_c = storage->get_data<event_cluster>("ImageClusterHit");
    auto ev_hit_02 = storage->get_data<event_hit>("hit02");

    if ( !ev_hit || !ev_hit->size()  ){
      std::cout<<"No hits here..."<<std::endl;
      return false ;
      } 

    if ( !ev_mcs || !ev_mcs->size()  ){
      std::cout<<"No mcshower here..."<<std::endl;
      return false ;
      } 

    if ( !ev_s || !ev_s->size()  ){
      std::cout<<"No shower here..."<<std::endl;
      return false ;
      } 

    if ( !ev_c || !ev_c->size()  ){
      std::cout<<"No cls here..."<<std::endl;
      return false ;
      } 

    int max_i = 0;
    int max_hits =0;
    auto ev_ass_c = storage->get_data<larlite::event_ass>("ImageClusterHit");
    auto const& ass_hit_v = ev_ass_c->association(ev_c->id(), ev_hit_02->id());


    for (size_t i = 0; i< ev_c->size(); ++i) {
      auto c = ev_c->at(i);
      if ( c.NHits() > max_hits ){
        max_hits = c.NHits();
	max_i = i;
           }
       }

     _sum_reco_int = 0;

     for (size_t i = 0; i< ass_hit_v.at(max_i).size(); ++i) {
         auto hit_id = ass_hit_v.at(max_i).at(i);
	 auto ihit = ev_hit_02->at(hit_id);
         _sum_reco_int += ihit.Integral();
       }



     


    _reco_e = 0;
    _sum = 0;
    _sum_adj = 0;
    for (size_t hindex = 0; hindex < ev_hit->size(); ++hindex) {
    
      auto const& h = (*ev_hit)[hindex];
      if ( h.WireID().Plane != 2) continue ;
      auto clocktick = larutil::DetectorProperties::GetME()->SamplingRate() * 1.e-3; //time sample in usec
      float lifetime_corr = exp( h.PeakTime() * clocktick / 1e6);
      float electrons = h.Integral() * 200 ; //mcc8 value

      float dQ = electrons * lifetime_corr * 23.6 * 1e-6 ; 
      float dE = dQ / 0.577 ; // 0.62 -> recomb factor

      _sum += h.Integral() ;

      _sum_adj += dE;

    }// for all h


  _mc_e = ev_mcs->at(0).Start().E() ;

  for (auto const & s : *ev_s ){
    if ( s.Energy(2) > _reco_e )
      _reco_e = s.Energy(2) ;
    }

  _tree->Fill();

  std::cout<<"MCShower energy : " <<ev_mcs->at(0).Start().E()<<std::endl ;
  std::cout<<"Reco energy : "<<_reco_e<<std::endl;
  
    return true;
  }

  bool GammaEnergy::finalize() {
    
    _tree->Write();
  
    return true;
  }

}
#endif
