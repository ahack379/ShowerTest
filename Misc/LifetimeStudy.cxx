#ifndef LARLITE_LIFETIMESTUDY_CXX
#define LARLITE_LIFETIMESTUDY_CXX

#include "LifetimeStudy.h"
#include "DataFormat/track.h"
#include "DataFormat/shower.h"
#include "DataFormat/cluster.h"
#include "DataFormat/opflash.h"
#include "DataFormat/hit.h"
#include "DataFormat/event_ass.h"

#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"

namespace larlite {

  bool LifetimeStudy::initialize() {

    if ( !_trk_tree ){
      _trk_tree = new TTree("trk_tree","trk_tree");
      _trk_tree->Branch("flash_t",&_flash_t,"flash_t/F");
      _trk_tree->Branch("flash_pe",&_flash_pe,"flash_pe/F");
      _trk_tree->Branch("hit_t",&_hit_t,"hit_t/F");
      _trk_tree->Branch("hit_charge",&_hit_charge,"hit_charge/F");
    }

    if ( !_shr_tree ){
      _shr_tree = new TTree("shr_tree","shr_tree");
      _shr_tree->Branch("flash_t",&_flash_t,"flash_t/F");
      _shr_tree->Branch("flash_pe",&_flash_pe,"flash_pe/F");
      _shr_tree->Branch("hit_shr_t",&_hit_shr_t,"hit_shr_t/F");
      _shr_tree->Branch("hit_shr_charge",&_hit_shr_charge,"hit_shr_charge/F");
    }

    if ( !_all_tree ){
      _all_tree = new TTree("all_tree","all_tree");
      _all_tree->Branch("flash_t",&_flash_t,"flash_t/F");
      _all_tree->Branch("flash_pe",&_flash_pe,"flash_pe/F");
      _all_tree->Branch("hit_t",&_hit_t,"hit_t/F");
      _all_tree->Branch("hit_charge",&_hit_charge,"hit_charge/F");
      _all_tree->Branch("is_shower",&_is_shower,"is_shower/I");
      _all_tree->Branch("pdg",&_pdg,"pdg/I");
      _all_tree->Branch("mom_pdg",&_mom_pdg,"mom_pdg/I");
    }

    return true;
  }
  
  bool LifetimeStudy::analyze(storage_manager* storage) {

    auto ev_flash = storage->get_data<larlite::event_opflash>("simpleFlashBeam");
    if ( ev_flash->size() ){
      int max_it = -1;
      int max_pe = -1;
      for( int ff = 0; ff < ev_flash->size(); ff++){
        auto f = ev_flash->at(ff);
        if ( f.TotalPE() > max_pe ) {
          max_pe = f.TotalPE();
          max_it = ff;
        }
      }
     if ( max_it != -1 ){
      auto f = ev_flash->at(max_it);
      _flash_t = f.Time();
      _flash_pe = f.TotalPE();
     }
    }

    // Now get Mccluster info
    auto ev_ass = storage->get_data<larlite::event_ass>("mccluster");
    auto const& ass_keys = ev_ass->association_keys();

    if ( ass_keys.size() == 0 ) return false;

    larlite::event_cluster *ev_mcclus = nullptr;
    auto ass_hit_clus_v = storage->find_one_ass( ass_keys[0].first, ev_mcclus, ev_ass->name() );

    larlite::event_hit *ev_mchit = nullptr;
    auto ass_mcclus_v = storage->find_one_ass( ass_keys[0].second, ev_mchit, ev_ass->name() );

    if (ass_hit_clus_v.size() == 0){
      std::cout << "No hit ass! exit" << std::endl;
      return false;
    }
    if (ass_mcclus_v.size() == 0){
      std::cout << "No mcclus ass! exit" << std::endl;
      return false;
    }

    auto ev_mct = storage->get_data<event_mctrack>("mcreco");
    if ( !ev_mct || !ev_mct->size() ) {std::cout<<"No MCTrack!" <<std::endl ; return false; }

    auto ev_mcs = storage->get_data<event_mcshower>("mcreco") ;
    if ( !ev_mcs || !ev_mcs->size() ) {std::cout<<"No MCShower!" <<std::endl ; return false; }

    auto ev_mcc = storage->get_data<event_cluster>("mccluster");

    auto ev_hit = storage->get_data<larlite::event_hit>("gaushit");

    // Fill map with hits from mccluster : clusterID
    for( int i = 0; i < ass_mcclus_v.size(); i ++ ){

      auto cid = ev_mcc->at(i) ;
      if ( cid.View() != 2 ) continue;

      _is_shower = cid.StartOpeningAngle() ; // opening angle set to track (0) or shower(1) in mccluster builder
      auto ts_index = cid.Width() ;           // width set to carry mct/s index

      if( !_is_shower ){
         auto mct = ev_mct->at(ts_index) ;
         if ( mct.Origin() != 1 ) continue;
         _pdg = mct.PdgCode();
         _mom_pdg = mct.MotherPdgCode();
       }
       else{
         auto mcs = ev_mcs->at(ts_index) ;
         if ( mcs.Origin() != 1 ) continue;
         _pdg = mcs.PdgCode();
         _mom_pdg = mcs.MotherPdgCode();
       }

      for ( int j = 0; j < ass_mcclus_v[i].size(); j++ ){

        auto hid = ass_mcclus_v[i][j];
        auto h = ev_hit->at(hid);

        _hit_t = h.PeakTime(); 
        _hit_charge = h.Integral();
        _all_tree->Fill();
      }    
    }    

    auto ev_tagged_trk = storage->get_data<event_track>("numuCC_track");
    if ( !ev_tagged_trk || !ev_tagged_trk->size() ){ std::cout<<"No Tagged Track!" <<std::endl ; return false; }

    auto tagged_trk = ev_tagged_trk->at(0) ;

    // Fill multiplicity info + find ID of pandora track that is numuCC_track
    auto ev_trk = storage->get_data<event_track>("pandoraNu");
    int min_trk_dist = 1e9;
    int min_trk_dist_it = -1;

    for ( int ii = 0; ii < ev_trk->size(); ii++){

        auto t = ev_trk->at(ii);
        auto st = t.Vertex() ;
        auto end = t.End() ;
        auto tag_end = tagged_trk.End() ;

        auto dist = sqrt( pow(tag_end.X() - end.X(),2) + pow(tag_end.Y() - end.Y(),2) + pow(tag_end.Z() - end.Z(),2) );
        if ( dist < min_trk_dist ){
          min_trk_dist = dist ;
          min_trk_dist_it = ii ;
        }
    }


    if (!ev_hit || ev_hit->size() == 0){
      std::cout << "No hits! exit" << std::endl;
      return false;
    }   

    auto ev_hit_cosRem = storage->get_data<event_hit>("pandoraCosmicHitRemoval"); 

    if ( !ev_hit_cosRem || ev_hit_cosRem->size() == 0 ) {
      std::cout << "No such hits associated to track! " << std::endl;
      return false;
    }

    auto ev_hr_ass = storage->get_data<larlite::event_ass>("pandoraNu"); 

    if ( !ev_hr_ass || ev_hr_ass->size() == 0 ) {
      std::cout << "No such association! " << std::endl;
      return false;
    }

    auto const& ass_hit_v = ev_hr_ass->association(ev_trk->id(), ev_hit_cosRem->id());

    if ( ass_hit_v.size() == 0) {
      std::cout << "No ass from track => hit! " << std::endl;
      return false;
    }

    // This particular headache is necessary because there are no gaushit associations
    // to pandoraNu tracks, only pandoraNuCosmicRemoval associations. 
     std::vector<int> tag_trk_gaushit_v;
     for(int i = 0; i < ass_hit_v.at(min_trk_dist_it).size(); i++){

        _hit_t = -999 ;
        _hit_charge = -999; 

        auto hid = ass_hit_v.at(min_trk_dist_it).at(i) ;
        auto h = ev_hit_cosRem->at(hid);
        if ( h.WireID().Plane != 2 ) continue;

        for(int j = 0; j < ev_hit->size(); j++){
          auto hj = ev_hit->at(j) ; 
          if ( hj.WireID().Plane != 2 ) continue;

          if ( hj.PeakTime() == h.PeakTime() && hj.WireID().Wire == h.WireID().Wire ){
            tag_trk_gaushit_v.emplace_back(j);
            _hit_t = hj.PeakTime() ;
            _hit_charge = hj.Integral() ;
            _trk_tree->Fill();
          }
        }
     }

     // Got tracks now get showers
     auto const & ev_clus = storage->get_data<event_cluster>("ImageClusterHit");
     auto const & ev_ass_c = storage->get_data<larlite::event_ass>("ImageClusterHit");
     auto const & ass_imageclus_v = ev_ass_c->association(ev_clus->id(), ev_hit->id());

     auto ev_s = storage->get_data<event_shower>("pi0_candidate_showers");
     if( !ev_s || !ev_s->size() ){ std::cout<<"Not enough reco'd showers..." <<std::endl; return false; }    

     // Get the association from shower -> cluster
     auto ev_ass_s = storage->get_data<larlite::event_ass>("pi0_candidate_showers");
     auto const& ass_showerreco_v = ev_ass_s->association(ev_s->id(), ev_clus->id());
     
     // Loop over showers associated to pi0 candidates
     for (size_t i = 0; i < ass_showerreco_v.size(); i++ ){
       for (size_t j = 0; j < ass_showerreco_v.at(i).size(); j++ ){

         auto clus_id = ass_showerreco_v.at(i).at(j); 
         auto iclus = ev_clus->at(clus_id);
     
         int plane = iclus.Plane().Plane ;
         if ( plane != 2 ) continue;

         // Loop through all hits associared to the cluster 
         for ( int k = 0; k < ass_imageclus_v.at(clus_id).size(); k++ ){

           auto hid = ass_imageclus_v.at(clus_id).at(k) ;
           auto h = ev_hit->at(hid);
           _hit_shr_t = h.PeakTime() ;
           _hit_shr_charge = h.Integral() ;
           _shr_tree->Fill(); 
         }
       }
     }

  
    return true;
  }

  bool LifetimeStudy::finalize() {

    if(_fout) { _fout->cd(); _trk_tree->Write(); _shr_tree->Write(); _all_tree->Write() ;}
  
    return true;
  }

}
#endif
