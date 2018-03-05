#ifndef LARLITE_POSTSEL2SHOWERQUALITY_CXX
#define LARLITE_POSTSEL2SHOWERQUALITY_CXX

#include "ShowerRecoEff.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/track.h"
#include "DataFormat/shower.h"
#include "DataFormat/vertex.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

#include "LArUtil/GeometryHelper.h"

namespace larlite {

  bool ShowerRecoEff::initialize() {    

    _event = -1; 

    _SCE = new larutil::SpaceChargeMicroBooNE();
    _time2cm = larutil::GeometryHelper::GetME()->TimeToCm();

    if ( !_tree ){
      _tree = new TTree("tree","");
      _tree->Branch("event",&_event,"event/I");

      _tree->Branch("purity",&_purity,"purity/F");
      _tree->Branch("complete",&_complete,"complete/F");
      _tree->Branch("cw_purity",&_cw_purity,"cw_purity/F");
      _tree->Branch("cw_complete",&_cw_complete,"cw_complete/F");
      _tree->Branch("origin",&_origin,"origin/F");
      _tree->Branch("type",&_type,"type/F");
      _tree->Branch("from_pi0",&_from_pi0,"from_pi0/B");

      _tree->Branch("reco_e",&_reco_e,"reco_e/F");
      _tree->Branch("st_x",&_st_x,"st_x/F");
      _tree->Branch("st_y",&_st_y,"st_y/F");
      _tree->Branch("st_z",&_st_z,"st_z/F");
      _tree->Branch("dir_x",&_dir_x,"dir_x/F");
      _tree->Branch("dir_y",&_dir_y,"dir_y/F");
      _tree->Branch("dir_z",&_dir_z,"dir_z/F");


      _tree->Branch("mc_e",&_mc_e,"mc_e/F");
      _tree->Branch("mc_detProf_e",&_mc_detProf_e,"_mc_detProf_e/F");
      _tree->Branch("mc_st_x",&_mc_st_x,"mc_st_x/F");
      _tree->Branch("mc_st_y",&_mc_st_y,"mc_st_y/F");
      _tree->Branch("mc_st_z",&_mc_st_z,"mc_st_z/F");
      _tree->Branch("mc_dir_x",&_mc_dir_x,"mc_dir_x/F");
      _tree->Branch("mc_dir_y",&_mc_dir_y,"mc_dir_y/F");
      _tree->Branch("mc_dir_z",&_mc_dir_z,"mc_dir_z/F");
      _tree->Branch("mc_detProf_x",&_mc_detProf_x,"mc_detProf_x/F");
      _tree->Branch("mc_detProf_y",&_mc_detProf_y,"mc_detProf_y/F");
      _tree->Branch("mc_detProf_z",&_mc_detProf_z,"mc_detProf_z/F");
      _tree->Branch("mc_dir_sce_corr_x",&_mc_dir_sce_corr_x,"mc_dir_sce_corr_x/F");
      _tree->Branch("mc_dir_sce_corr_y",&_mc_dir_sce_corr_y,"mc_dir_sce_corr_y/F");
      _tree->Branch("mc_dir_sce_corr_z",&_mc_dir_sce_corr_z,"mc_dir_sce_corr_z/F");
      _tree->Branch("mc_detProf_sce_corr_x",&_mc_detProf_sce_corr_x,"mc_detProf_sce_corr_x/F");
      _tree->Branch("mc_detProf_sce_corr_y",&_mc_detProf_sce_corr_y,"mc_detProf_sce_corr_y/F");
      _tree->Branch("mc_detProf_sce_corr_z",&_mc_detProf_sce_corr_z,"mc_detProf_sce_corr_z/F");

      _tree->Branch("n_true_showers",&_n_true_showers,"n_true_showers/F");
      _tree->Branch("n_recod_true_showers",&_n_recod_true_showers,"n_recod_true_showers/F");
      _tree->Branch("low_shr_e",&_low_shr_e,"low_shr_e/F");
      _tree->Branch("high_shr_e",&_high_shr_e,"high_shr_e/F");

   }
   _out_of_av = 0;
   _tot_shr = 0;

    return true;
  }

  void ShowerRecoEff::clear(){

    _purity = 0.;
    _complete = 0.;
    _cw_purity = 0.;
    _cw_complete = 0.;

    _mc_e = -999;
    _mc_detProf_e = -999;
    _origin = -1;
    _type = -1 ; // 1 is shower, 0 is track
    _from_pi0= false ; 

    _st_x  = -999;
    _st_y  = -999;
    _st_z  = -999;
    _dir_x  = -999;
    _dir_y  = -999;
    _dir_z  = -999;

    _mc_st_x = -999;
    _mc_st_y = -999;
    _mc_st_z = -999;
    _mc_dir_x = -999;
    _mc_dir_y = -999;
    _mc_dir_z = -999;
    _mc_detProf_x = -999;
    _mc_detProf_y = -999;
    _mc_detProf_z = -999;
    _mc_dir_sce_corr_x = -999;
    _mc_dir_sce_corr_y = -999;
    _mc_dir_sce_corr_z = -999;
    _mc_detProf_sce_corr_x = -999;
    _mc_detProf_sce_corr_y = -999;
    _mc_detProf_sce_corr_z = -999;
    
    _n_true_showers = 0;
    _n_recod_true_showers = 0;

    _low_shr_e = 1e9; 
    _high_shr_e = -999;

  }
  
  bool ShowerRecoEff::analyze(storage_manager* storage) {

    _event++ ;
    std::cout<<"\n\nEVENT IS: "<<_event<<std::endl;
    clear();

    auto ev_hit = storage->get_data<larlite::event_hit>("gaushit");

    if (!ev_hit || ev_hit->size() == 0){
      std::cout << "No hits! exit" << std::endl;
      return false;
    }   

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
    }

    // Now get Mccluster info
    auto ev_mcc = storage->get_data<event_cluster>("mccluster");
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

    // Keep track of the charge-weighted hit count
    std::map<int,float> tot_mc_cw_hits_v ; 

    _mc_hit_map.clear();

    // Fill map with hits from mccluster : clusterID
    for( int i = 0; i < ass_mcclus_v.size(); i ++ ){

      auto cid = ev_mcc->at(i) ;
      if ( cid.View() != 2 ) continue;

      for ( int j = 0; j < ass_mcclus_v[i].size(); j++ ){

        auto hid = ass_mcclus_v[i][j];
        _mc_hit_map[hid] = i ; 

        auto h = ev_hit->at(hid);

        if ( tot_mc_cw_hits_v.find(i) == tot_mc_cw_hits_v.end() )
          tot_mc_cw_hits_v[i] = h.Integral() ;
        else
          tot_mc_cw_hits_v[i] += h.Integral() ;
      }
    }

    std::vector<int> pur_ctr_v ;
    std::vector<float> cw_pur_ctr_v ;

    int max_hits = -1;
    int max_cw_hits = -1;
    int max_cid = -1 ;
    float tot_reco_cw_hits = 0;

    // Get the association from cluster -> hit
    auto const & ev_clus = storage->get_data<event_cluster>("ImageClusterHit");
    auto const & ev_ass_c = storage->get_data<larlite::event_ass>("ImageClusterHit");
    auto const & ass_imageclus_v = ev_ass_c->association(ev_clus->id(), ev_hit->id());

    auto ev_mcs = storage->get_data<event_mcshower>("mcreco") ;
    if ( !ev_mcs || !ev_mcs->size() ) {std::cout<<"No MCShower!" <<std::endl ; return false; }

    auto ev_s = storage->get_data<event_shower>("showerreco");

    //std::cout<<"Reconstructed showers: "<<ev_s->size()<<std::endl;
    //if( !ev_s || !ev_s->size() ){ 
    //  std::cout<<"Not enough reco'd showers..." <<std::endl;
    //  return false;
    //}   

    auto temp_n_true_showers = 0;
    auto temp_n_recod_true_showers = 0;
    std::vector<int> all_mcs_v ;  
    std::vector<int> used_mcs_v ;  

    for ( int i = 0; i < ev_mcs->size(); i++){

      auto s = ev_mcs->at(i);
    
      // Need this line for neutrino interactions; don't for single particle
      if ( s.Origin() != 1 || s.MotherPdgCode() != 111 || s.AncestorPdgCode() != 111) continue;
      //if ( s.MotherPdgCode() != 111 || s.AncestorPdgCode() != 111) continue;
      
      auto x = s.Start().X() ;
      auto y = s.Start().Y() ;
      auto z = s.Start().Z() ;

      if ( x < 20 || x > 236.35 || y < -96.5 || y > 96.5 || z < 10 || z > 1026.7 ) {std::cout<<"OUT OF FV ! "<<std::endl; return false; }

      all_mcs_v.emplace_back(i);
      
      if ( s.DetProfile().E() < _low_shr_e )
        _low_shr_e = s.DetProfile().E() ;
      if ( s.DetProfile().E() > _high_shr_e )
        _high_shr_e = s.DetProfile().E() ;

      //std::cout<<"Mother: "<<s.AncestorPdgCode()<<std::endl ;

      temp_n_true_showers++;
      auto end = s.End() ;
      _tot_shr ++; 

      if( end.X() < 0 || end.X() > 256.35 || end.Y() > 116.5 || end.Y() < -116.5 || end.Z() < 0 || end.Z() > 1036.8 )
        _out_of_av ++;
    
    }

    if ( ev_s->size() ){   

      // Get the association from shower -> cluster
      auto ev_ass_s = storage->get_data<larlite::event_ass>("showerreco");
      auto const& ass_showerreco_v = ev_ass_s->association(ev_s->id(), ev_clus->id());
      // Loop over showers
      for (size_t i = 0; i < ass_showerreco_v.size(); i++ ){

        float mc_clus_e = 0.;
            int closest_mcs_id = -1 ;

        // Loop over clusters associated to this shower
        for (size_t j = 0; j < ass_showerreco_v.at(i).size(); j++ ){

              auto clus_id = ass_showerreco_v.at(i).at(j); 
              auto iclus = ev_clus->at(clus_id);
          
              int plane = iclus.Plane().Plane ;
              if ( plane != 2 ) continue;

              pur_ctr_v.clear();
              cw_pur_ctr_v.clear();
              pur_ctr_v.resize(ass_mcclus_v.size()+1,0) ;
              cw_pur_ctr_v.resize(ass_mcclus_v.size()+1,0) ;

              max_hits = -1;
              max_cw_hits = -1;
              max_cid = -1 ;
              tot_reco_cw_hits = 0;

              // Loop through all hits associared to the cluster 
              for ( int k = 0; k < ass_imageclus_v.at(clus_id).size(); k++ ){

                auto hid = ass_imageclus_v.at(clus_id).at(k) ; 
                auto h = ev_hit->at(hid);
                tot_reco_cw_hits += h.Integral() ;
                
                if ( _mc_hit_map.find(hid) != _mc_hit_map.end() ){

                  auto mcclus_id = _mc_hit_map[hid] ;

                  pur_ctr_v[mcclus_id]++ ; 
                  cw_pur_ctr_v[mcclus_id] += h.Integral() ; 

                  if( pur_ctr_v[ mcclus_id] > max_hits ){
                    max_hits = pur_ctr_v[mcclus_id];
                    max_cid = mcclus_id ; 
                    max_cw_hits = cw_pur_ctr_v[mcclus_id] ;
                  }
                }
                    else {
                      auto mcclus_id = ass_mcclus_v.size() ;
                  pur_ctr_v[mcclus_id]++ ; 
                  cw_pur_ctr_v[mcclus_id] += h.Integral() ; 
                  if( pur_ctr_v[mcclus_id] > max_hits ){
                    max_hits = pur_ctr_v[mcclus_id];
                    max_cid = mcclus_id ; 
                    max_cw_hits = cw_pur_ctr_v[mcclus_id] ;
                  }
                    }
              }

              if ( max_cid != ass_mcclus_v.size() && max_cid != -1 ){

                auto mcclus = ev_mcc->at(max_cid) ;
                auto clocktick = larutil::DetectorProperties::GetME()->SamplingRate() * 1.e-3; //time sample in usec
        	
                // Store true shower detprofile energy
                for ( auto const & mcc_hid : ass_mcclus_v[max_cid] ){
                  auto mch = ev_hit->at(mcc_hid) ;
                  float lifetime_corr = exp( mch.PeakTime() * clocktick / 1.e20);
                  float electrons = mch.Integral() * 198.; //mcc8 value
                  float dQ = electrons * lifetime_corr * 23.6 * 1e-6 ;
                  float dE = dQ / 0.577 ; // 0.62 -> recomb factor
                  mc_clus_e += dE ;
                }

                float closest_e = 1e9 ;

                // Find mcs this cluster belongs to in order to store the true shower energy
                for ( int i = 0; i < ev_mcs->size(); i++ ) { 

        	        auto s = ev_mcs->at(i) ;
        	        if(s.Origin() != 1 || s.MotherPdgCode() != 111 ) continue;
        	        //if( s.MotherPdgCode() != 111 ) continue;
                        
                         auto e = fabs(mc_clus_e - s.DetProfile().E()) ;

        	        if ( e < closest_e ){
        	          closest_e  = e;
        	          closest_mcs_id = i ;
        	        }
                }

                used_mcs_v.emplace_back(closest_mcs_id);

                auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
                auto tot_reco_hits = ass_imageclus_v[clus_id].size();
                
                _purity   = float(max_hits) / tot_reco_hits ;
                _complete = float(max_hits) / tot_mc_hits ;

                _cw_purity   = float(max_cw_hits) / tot_reco_cw_hits ;
                _cw_complete = float(max_cw_hits) / tot_mc_cw_hits_v[max_cid]; 

                //auto mcclus = ev_mcc->at(max_cid) ;
                _origin = mcclus.Width() ; 
                _type   = mcclus.StartOpeningAngle() ; // Recall I've set this to track (0) or shower(1) in mccluster builder
                _from_pi0 = mcclus.IsMergedCluster() ; 
             }
           }

           auto ishr = ev_s->at(i);

           _reco_e = ishr.Energy(2);
           _st_x = ishr.ShowerStart().X();
           _st_y = ishr.ShowerStart().Y();
           _st_z = ishr.ShowerStart().Z();
           _dir_x = ishr.Direction().X();
           _dir_y = ishr.Direction().Y();
           _dir_z = ishr.Direction().Z();

           if ( closest_mcs_id != -1 )
             _mc_detProf_e = ev_mcs->at(closest_mcs_id).DetProfile().E() ; 
             //_mc_detProf_e =  mc_clus_e;

           if ( _origin == 1 && _type == 1 && _from_pi0 ){ // This one is for events with neutrino interactions 
           //if ( _from_pi0 ){ // This one is for single particle sample
             
             _n_recod_true_showers = 1;
             _n_true_showers = 1 ;
             temp_n_recod_true_showers ++; 

             _tree->Fill();    
           }
        }
     }

     std::cout<<"reco and true : "<<temp_n_recod_true_showers<<", "<<temp_n_true_showers<<std::endl ;
     if ( temp_n_recod_true_showers < temp_n_true_showers ){

       for ( int i = 0 ; i < all_mcs_v.size() ; i++ ){

         int id = all_mcs_v.at(i);
	 if( std::find(used_mcs_v.begin(),used_mcs_v.end(),id) != used_mcs_v.end() )
	   continue;
	 else{

	   _n_recod_true_showers = 0;
	   _n_true_showers = 1;
	   _mc_detProf_e = ev_mcs->at(id).DetProfile().E() ;

	   std::cout<<"Filling MCCC "<<std::endl ;
	   _tree->Fill();
	 
	 }
       
       }
     }
    
    return true;
  }

  bool ShowerRecoEff::finalize() {

    std::cout<<"Out of AV: "<<_out_of_av<<", "<<_tot_shr<<std::endl ;

    if ( _fout ){
      _fout->cd();
      _tree->Write();
    }
  
    return true;
  }

}
#endif
