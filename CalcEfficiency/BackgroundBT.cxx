#ifndef LARLITE_BACKGROUNDBT_CXX
#define LARLITE_BACKGROUNDBT_CXX

#include "BackgroundBT.h"
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

  bool BackgroundBT::initialize() {    

    _event = -1; 

    _n_noise = 0;    // 0
    _n_cosmic = 0;   // 1         
    _n_cc1pi0 = 0;   // 2
    _n_nc1pi0 = 0;   // 3
    _n_cc1pi0_outFV = 0;  // 4
    _n_multpi0 = 0;  // 5
    _n_nue = 0;      // 6 
    _n_antimu  = 0;   // 7
    _n_cccex = 0;    // 8
    _n_nccex = 0;    // 9
    _n_ccother = 0;  // 12
    _n_ncother = 0;  // 13
    _n_other = 0;    // 14 

    // After thoughts 
    _n_gamma = 0 ; // 10
    _n_kaondecay = 0 ; // 11

    _bkgd_v = { "Other","Cosmic","CC1pi0","NC1pi0","CC1pi0_outFV","Multpi0", "nue","antinumu","cccex","nccex","ccgamma","kaondecay","ccother","ncother" };

    if ( !_tree ){
      _tree = new TTree("tree","");
      _tree->Branch("event",&_event,"event/I");
      _tree->Branch("bkgd_id",&_bkgd_id,"bkgd_id/I");
      _tree->Branch("nu_mode",&_nu_mode,"nu_mode/I");
      _tree->Branch("nshrs",&_nshrs,"nshrs/I");

      _tree->Branch("vtx_x",&_vtx_x,"vtx_x/F");
      _tree->Branch("vtx_y",&_vtx_y,"vtx_y/F");
      _tree->Branch("vtx_z",&_vtx_z,"vtx_z/F");
      _tree->Branch("mc_vtx_x",&_mc_vtx_x,"mc_vtx_x/F");
      _tree->Branch("mc_vtx_y",&_mc_vtx_y,"mc_vtx_y/F");
      _tree->Branch("mc_vtx_z",&_mc_vtx_z,"mc_vtx_z/F");

      _tree->Branch("mu_angle",&_mu_angle,"mu_angle/F");
      _tree->Branch("mu_len",&_mu_len,"mu_len/F");
      _tree->Branch("mu_startx",&_mu_startx,"mu_startx/F");
      _tree->Branch("mu_starty",&_mu_starty,"mu_starty/F");
      _tree->Branch("mu_startz",&_mu_startz,"mu_startz/F");
      _tree->Branch("mu_endx",&_mu_endx,"mu_endx/F");
      _tree->Branch("mu_endy",&_mu_endy,"mu_endy/F");
      _tree->Branch("mu_endz",&_mu_endz,"mu_endz/F");
      _tree->Branch("mu_phi",&_mu_phi,"mu_phi/F");
      _tree->Branch("mu_mom",&_mu_mom,"mu_mom/F");
      _tree->Branch("mult",&_mult,"mult/F");
      _tree->Branch("mu_purity",&_mu_purity,"mu_purity/F");
      _tree->Branch("mu_complete",&_mu_complete,"mu_complete/F");
      _tree->Branch("mu_cw_purity",&_mu_cw_purity,"mu_cw_purity/F");
      _tree->Branch("mu_cw_complete",&_mu_cw_complete,"mu_cw_complete/F");

      _tree->Branch("pi0_mass",&_pi0_mass,"pi0_mass/F");
      _tree->Branch("pi0_oangle",&_pi0_oangle,"pi0_oangle/F");
      _tree->Branch("pi0_IP",&_pi0_IP,"pi0_IP/F");
      _tree->Branch("pi0_mom",&_pi0_mom,"pi0_mom/F");
      _tree->Branch("pi0_low_shrE",&_pi0_low_shrE,"pi0_low_shrE/F");
      _tree->Branch("pi0_high_shrE",&_pi0_high_shrE,"pi0_high_shrE/F");
      _tree->Branch("pi0_low_radL",&_pi0_low_radL,"pi0_low_radL/F");
      _tree->Branch("pi0_high_radL",&_pi0_high_radL,"pi0_high_radL/F");

      _tree->Branch("pi0_low_purity",&_pi0_low_purity,"pi0_low_purity/F");
      _tree->Branch("pi0_high_purity",&_pi0_high_purity,"pi0_high_purity/F");
      _tree->Branch("pi0_low_complete",&_pi0_low_complete,"pi0_low_complete/F");
      _tree->Branch("pi0_high_complete",&_pi0_high_complete,"pi0_high_complete/F");

      _tree->Branch("pi0_low_cw_purity",&_pi0_low_cw_purity,"pi0_low_cw_purity/F");
      _tree->Branch("pi0_high_cw_purity",&_pi0_high_cw_purity,"pi0_high_cw_purity/F");
      _tree->Branch("pi0_low_cw_complete",&_pi0_low_cw_complete,"pi0_low_cw_complete/F");
      _tree->Branch("pi0_high_cw_complete",&_pi0_high_cw_complete,"pi0_high_cw_complete/F");


      _tree->Branch("gamma_E",&_gamma_E,"gamma_E/F");
      _tree->Branch("gamma_RL",&_gamma_RL,"gamma_RL/F");
      _tree->Branch("gamma_purity",&_gamma_purity,"gamma_purity/F");
      _tree->Branch("gamma_complete",&_gamma_complete,"gamma_complete/F");
      _tree->Branch("gamma_cw_purity",&_gamma_cw_purity,"gamma_cw_purity/F");
      _tree->Branch("gamma_cw_complete",&_gamma_cw_complete,"gamma_cw_complete/F");
   }

    return true;
  }

  void BackgroundBT::clear(){

    _bkgd_id = -1 ;
    _nu_mode = -1 ;
    _nshrs = -1 ;
    _mult    = 0;
    _vtx_x   = -999;
    _vtx_y   = -999;
    _vtx_z   = -999;
    _mc_vtx_x   = -999;
    _mc_vtx_y   = -999;
    _mc_vtx_z   = -999;

    _mu_angle = -999;
    _mu_len   = -999;
    _mu_startx = -999;
    _mu_starty = -999;
    _mu_startz = -999;
    _mu_endx = -999;
    _mu_endy = -999;
    _mu_endz = -999;
    _mu_phi  = -999;
    _mu_mom  = -999;
    _mu_complete  = 0.;
    _mu_purity  = 0.;
    _mu_cw_complete  = 0.;
    _mu_cw_purity  = 0.;

    _pi0_mass = -999;
    _pi0_oangle = -999;
    _pi0_IP = -999;
    _pi0_mom = -999;
    _pi0_low_shrE = -999;
    _pi0_high_shrE = -999;
    _pi0_low_radL = -999;
    _pi0_high_radL = -999;

    _pi0_high_purity = 0.;
    _pi0_high_complete = 0.;
    _pi0_low_purity = 0.;
    _pi0_low_complete = 0.;

    _pi0_high_cw_purity = 0.;
    _pi0_high_cw_complete = 0.;
    _pi0_low_cw_purity = 0.;
    _pi0_low_cw_complete = 0.;

    _gamma_E = -999;
    _gamma_RL = -999;
    _gamma_purity = 0.;
    _gamma_complete = 0.;
    _gamma_cw_purity = 0.;
    _gamma_cw_complete = 0.;

  }
  
  bool BackgroundBT::analyze(storage_manager* storage) {

    _event++ ;
    clear();

    auto ev_vtx= storage->get_data<event_vertex>("numuCC_vertex"); 
    if(!ev_vtx || !ev_vtx->size() ) {
      std::cout<<"Event has no recovertex info "<<std::endl;
      return false;
    }

    auto vtx = ev_vtx->at(0); 
    _vtx_x = vtx.X();
    _vtx_y = vtx.Y();
    _vtx_z = vtx.Z();

    auto ev_tagged_trk = storage->get_data<event_track>("numuCC_track");
    if ( !ev_tagged_trk || !ev_tagged_trk->size() ){ std::cout<<"No Tagged Track!" <<std::endl ; return false; }

    auto tagged_trk = ev_tagged_trk->at(0) ;

    // Fill track information
    _mu_phi = tagged_trk.Phi();
    _mu_angle = cos(tagged_trk.Theta());
    _mu_len =   tagged_trk.Length(0); // Calculates the length from point 0 to end
    _mu_startx = tagged_trk.Vertex().X(); 
    _mu_starty = tagged_trk.Vertex().Y(); 
    _mu_startz = tagged_trk.Vertex().Z(); 
    _mu_endx = tagged_trk.End().X(); 
    _mu_endy = tagged_trk.End().Y(); 
    _mu_endz = tagged_trk.End().Z(); 

    _mu_mom  = tagged_trk.VertexMomentum() ;

    // Adjust for pandora bug
    std::vector<double> dir = { (_mu_endx - _mu_startx) / _mu_len,
                                (_mu_endy - _mu_starty) / _mu_len,
                                (_mu_endz - _mu_startz) / _mu_len };

    auto dir_start = tagged_trk.VertexDirection();
    std::vector<double> other_dir = { dir_start.X(), dir_start.Y(), dir_start.Z() };  

    float dotProd = dir.at(0) * other_dir.at(0) + dir.at(1) * other_dir.at(1) +  dir.at(2) * other_dir.at(2) ;

    if( dotProd < 0 ) { 
       TVector3 new_dir(-dir_start.X(),-dir_start.Y(),-dir_start.Z());
       _mu_angle = cos(new_dir.Theta());
       _mu_phi = new_dir.Phi();
    }   

    // Fill multiplicity info + find ID of pandora track that is numuCC_track
    auto ev_trk = storage->get_data<event_track>("pandoraNu");
    int min_trk_dist = 1e9;
    int min_trk_dist_it = -1;

    //for ( auto const & t : *ev_trk ){
    for ( int ii = 0; ii < ev_trk->size(); ii++){
        auto t = ev_trk->at(ii);
        auto st = t.Vertex() ;
        auto end = t.End() ;

        auto dist_st = sqrt( pow(st.X() - vtx.X(),2) + pow(st.Y() - vtx.Y(),2) + pow(st.Z() - vtx.Z(),2) );
        auto dist_end = sqrt( pow(end.X() - vtx.X(),2) + pow(end.Y() - vtx.Y(),2) + pow(end.Z() - vtx.Z(),2) );

        if (dist_st < 3 || dist_end < 3)
          _mult ++ ;

        auto tag_st = tagged_trk.Vertex() ;
        auto dist = sqrt( pow(tag_st.X() - vtx.X(),2) + pow(tag_st.Y() - vtx.Y(),2) + pow(tag_st.Z() - vtx.Z(),2) );
        if ( dist < min_trk_dist ){
	  min_trk_dist = dist ;
	  min_trk_dist_it = ii ;
	}

    }

    auto ev_hit = storage->get_data<larlite::event_hit>("gaushit");

    if (!ev_hit || ev_hit->size() == 0){
      std::cout << "No hits! exit" << std::endl;
      return false;
    }   

    std::vector<int> pur_ctr_v ;
    std::vector<float> cw_pur_ctr_v ;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Want to be able to access the origin of the tagged muon. Thus, need to find it, and 
    // Ask for its origin.  Need to match to MCtrack to do this
    if ( _mc_sample ){

      auto ev_mctrk = storage->get_data<event_mctrack>("mcreco");
      if ( !ev_mctrk || !ev_mctrk->size() ) {std::cout<<"No MCTrack!" <<std::endl ; return false; }

      auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
      if(!ev_mctruth || !ev_mctruth->size() ) {
        std::cout<<"Event has no mctruth info "<<std::endl;
        return false;
      }

      auto & truth = ev_mctruth->at(0);
      auto & nu  = truth.GetNeutrino();

      _nu_mode = nu.Mode();
      if (_nu_mode == 10) _nu_mode = 4;

      auto traj = nu.Nu().Trajectory();
      _mc_vtx_x = traj.at(traj.size() - 1).X();
      _mc_vtx_y = traj.at(traj.size() - 1).Y();
      _mc_vtx_z = traj.at(traj.size() - 1).Z();

      auto vtx_diff = sqrt(pow(_mc_vtx_x - _vtx_x,2) + pow(_mc_vtx_y - _vtx_y,2) + pow(_mc_vtx_z - _vtx_z,2));

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


      // Keep track of the charge-weighted hit count
      std::map<int,float> tot_mc_cw_hits_v ; 

      _mc_hit_map.clear();

      // Fill map with hits from mccluster : clusterID
      for( int i = 0; i < ass_mcclus_v.size(); i ++ ){

        //std::cout<<"MCCluster id + hits : "<<i<<", "<<ass_mcclus_v[i].size()<<std::endl ;
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

      // Get association to trk => hit and hit => trk
      auto const& ass_hit_v = ev_hr_ass->association(ev_trk->id(), ev_hit_cosRem->id());

      if ( ass_hit_v.size() == 0) {
        std::cout << "No ass from track => hit! " << std::endl;
        return false;
      }

      // for each reco cluster, find the origin of all hits and calc purity/completeness 
      pur_ctr_v.resize(ass_mcclus_v.size(),0) ;
      cw_pur_ctr_v.resize(ass_mcclus_v.size(),0) ;

      int max_hits = -1;
      int max_cw_hits = -1;
      int max_cid = -1 ;
      float tot_reco_cw_hits = 0;
            
      for(int i = 0; i < ass_hit_v.at(min_trk_dist_it).size(); i++){
         auto hid = ass_hit_v.at(min_trk_dist_it).at(i) ;
         auto h = ev_hit->at(hid);
         if ( h.WireID().Plane != 2 ) continue;

         tot_reco_cw_hits += h.Integral() ;
         
         if ( _mc_hit_map.find(hid) != _mc_hit_map.end() ){

           auto mcclus_id = _mc_hit_map[hid] ;

           pur_ctr_v[mcclus_id]++ ; 
           cw_pur_ctr_v[mcclus_id] += h.Integral() ; 

           if( pur_ctr_v[mcclus_id] > max_hits ){
            
             max_hits = pur_ctr_v[mcclus_id];
             max_cid = mcclus_id ; 
             max_cw_hits = cw_pur_ctr_v[mcclus_id] ;
           }
         }
       }

       if ( max_cid != -1 ){

         auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
         auto tot_reco_hits = ass_hit_v.at(min_trk_dist_it).size();         
         _mu_purity   = float(max_hits) / tot_reco_hits ;
         _mu_complete = float(max_hits) / tot_mc_hits ;

         _mu_cw_purity   = float(max_cw_hits) / tot_reco_cw_hits ;
         _mu_cw_complete = float(max_cw_hits) / tot_mc_cw_hits_v[max_cid]; 
       }
      else {
         // Noise category
         _n_noise++;
         _bkgd_id = 0 ;
         _event_list.emplace_back(_event) ;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////
      /// Now count the number of backgrounds and signals
      ///////////////////////////////////////////////////////////////////////////////////////////////
      if( _bkgd_id == -1 ){

        auto ev_mcs = storage->get_data<event_mcshower>("mcreco") ;
        if ( !ev_mcs || !ev_mcs->size() ) {std::cout<<"No MCShower!" <<std::endl ; return false; }

        bool infv = true;
        if( _mc_vtx_x < 20 || _mc_vtx_x > 236.35 || _mc_vtx_y > 96.5 || _mc_vtx_y < -96.5 || _mc_vtx_z < 10 || _mc_vtx_z > 1026.8 )
          infv = false;

        auto parts = ev_mctruth->at(0).GetParticles();
        int n_pi0 = 0;
        int n_gamma = 0;

        for ( auto const & p : parts ){
          if( p.StatusCode() == 1 && p.PdgCode() == 111 )
            n_pi0++;

          if( p.StatusCode() == 1 && p.PdgCode() == 22 )
            n_gamma++;
        }   

       auto ev_mcc = storage->get_data<event_cluster>("mccluster");
       auto mcclus = ev_mcc->at(max_cid) ;
       // Remember that I've repurposed the cluster width variable to store info 
       // about whether this mccluster (which we've identified as the match to 
       // our reco particle at this stage) is neutrino or cosmic in origin 

        if( mcclus.Width() == 2){
          _n_cosmic++;
          _bkgd_id = 1; 
        }
        else if( abs(nu.Nu().PdgCode()) == 12 ){
          _n_nue ++ ; 
          _bkgd_id = 6 ;
        }
        else if( nu.Nu().PdgCode() == -14 ){
          _n_antimu ++ ; 
          _bkgd_id = 7 ;
        }
        // frmo here we can assume we have a muon neutrino
        else if( nu.CCNC() == 0 && n_pi0 == 1 && infv ) {
          _bkgd_id = 2;
          _n_cc1pi0++; 
        }
        else if( nu.CCNC() == 1 && n_pi0 > 0 ) {
          _bkgd_id = 3;
          _n_nc1pi0 ++; 
        }
        else if ( nu.CCNC() == 0 & n_pi0 == 1 && !infv ){
          _bkgd_id = 4 ;
          _n_cc1pi0_outFV ++; 
        }
        else if( nu.CCNC() == 0 && n_pi0 > 1 ) {
          _bkgd_id = 5;
          _n_multpi0 ++; 
        }
        else if( n_pi0 == 0 && n_gamma > 0 ){
          _bkgd_id = 10 ;
          _n_gamma++;
        }
        else{
          
          bool charge_ex = false; 
          bool kaon_decay = false; 
          for ( auto const & s : *ev_mcs ){
            if ( s.MotherPdgCode() == 111 && s.Origin() == 1 && abs(s.AncestorPdgCode()) == 211 ){
              charge_ex = true; 
              break;
            }

            if ( s.MotherPdgCode() == 111 && s.Origin() == 1 && abs(s.AncestorPdgCode()) == 321 ){
              kaon_decay = true; 
              break;
            }

          }

          if( charge_ex && nu.CCNC() == 0 ){
            _bkgd_id = 8;
            _n_cccex++;
          }
          else if( charge_ex && nu.CCNC() == 1 ) {
            _bkgd_id = 9;
            _n_nccex++;
          }
          else if( kaon_decay ){
            _bkgd_id = 11;
            _n_kaondecay++;
          }
          else if( !charge_ex && nu.CCNC() == 0 ){
            _bkgd_id = 12;
            _n_ccother++; 

          }
          else if( !charge_ex && nu.CCNC() == 1 ){
            _bkgd_id = 13;
            _n_ncother++; 
          }
          else {
            _bkgd_id = 14;
            _n_other ++;   
          }
        }
      }

      // Get the association from cluster -> hit
      auto const & ev_clus = storage->get_data<event_cluster>("ImageClusterHit");
      auto const & ev_ass_c = storage->get_data<larlite::event_ass>("ImageClusterHit");
      auto const & ass_imageclus_v = ev_ass_c->association(ev_clus->id(), ev_hit->id());

      if ( _get_pi0_info ){

        auto ev_s = storage->get_data<event_shower>("pi0_candidate_showers");
        if( !ev_s || !ev_s->size() ){ 
          std::cout<<"Not enough reco'd showers..." <<std::endl;
          return false;
         }   

        auto const& shr1 = ev_s->at(0) ;
        auto const& shr2 = ev_s->at(1);
        bool lowE_is_shr1 = shr1.Energy() < shr2.Energy() ? 1 : 0 ;

        // Get the association from shower -> cluster
        auto ev_ass_s = storage->get_data<larlite::event_ass>("pi0_candidate_showers");
        auto const& ass_showerreco_v = ev_ass_s->association(ev_s->id(), ev_clus->id());
        
        // Loop over showers
        for (size_t i = 0; i < ass_showerreco_v.size(); i++ ){

          float purity = 0., complete = 0., cw_purity = 0., cw_complete = 0.;

          // Loop over clusters associated to this shower
          for (size_t j = 0; j < ass_showerreco_v.at(i).size(); j++ ){

              auto clus_id = ass_showerreco_v.at(i).at(j); 
              auto iclus = ev_clus->at(clus_id);
          
              int plane = iclus.Plane().Plane ;
              if ( plane != 2 ) continue;

              pur_ctr_v.clear();
              cw_pur_ctr_v.clear();
              pur_ctr_v.resize(ass_mcclus_v.size(),0) ;
              cw_pur_ctr_v.resize(ass_mcclus_v.size(),0) ;

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
              }

              if ( max_cid != -1 ){

                auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
                auto tot_reco_hits = ass_imageclus_v[clus_id].size();
                
                purity   = float(max_hits) / tot_reco_hits ;
                complete = float(max_hits) / tot_mc_hits ;

                cw_purity   = float(max_cw_hits) / tot_reco_cw_hits ;
                cw_complete = float(max_cw_hits) / tot_mc_cw_hits_v[max_cid]; 
             }
           }
           if ( i == 0 ){
             if ( lowE_is_shr1 ){ 
               _pi0_low_purity = purity ;
               _pi0_low_complete = complete;
             }
             else{
               _pi0_high_purity = purity;
               _pi0_high_complete = purity;
             }
           }
           else if ( i == 1) {

             if ( lowE_is_shr1 ){ 
               _pi0_low_purity = purity ;
               _pi0_low_complete = complete;
             }
             else{
               _pi0_high_purity = purity;
               _pi0_high_complete = complete;
             }
           }
        }
      }

    if ( _get_single_shower_info ){
      auto ev_s = storage->get_data<event_shower>("pi0_1gamma_candidate_showers");
      if ( !ev_s ) return false;

      // Get the association from shower -> cluster
      auto ev_ass_s = storage->get_data<larlite::event_ass>("pi0_1gamma_candidate_showers");
      auto const& ass_showerreco_v = ev_ass_s->association(ev_s->id(), ev_clus->id());
      
      // Loop over clusters associated to this shower
      for (size_t j = 0; j < ass_showerreco_v.at(0).size(); j++ ){

        auto clus_id = ass_showerreco_v.at(0).at(j); 
        auto iclus = ev_clus->at(clus_id);
        
        int plane = iclus.Plane().Plane ;
        if ( plane != 2 ) continue;

        pur_ctr_v.clear();
        cw_pur_ctr_v.clear();
        pur_ctr_v.resize(ass_mcclus_v.size(),0) ;
        cw_pur_ctr_v.resize(ass_mcclus_v.size(),0) ;

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
        }

        if ( max_cid != -1 ){

          auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
          auto tot_reco_hits = ass_imageclus_v[clus_id].size();
          
          _gamma_purity   = float(max_hits) / tot_reco_hits ;
          _gamma_complete = float(max_hits) / tot_mc_hits ;

          _gamma_cw_purity   = float(max_cw_hits) / tot_reco_cw_hits ;
          _gamma_cw_complete = float(max_cw_hits) / tot_mc_cw_hits_v[max_cid]; 
        }
      }
    }
  }

    auto ev_shr = storage->get_data<event_shower>("showerreco");

    if ( ev_shr ) _nshrs = ev_shr->size(); 

    if ( _get_pi0_info ){

      auto ev_s = storage->get_data<event_shower>("pi0_candidate_showers");
      if( !ev_s || !ev_s->size() ){ //|| ev_s->size() < 2 ){
        std::cout<<"Not enough reco'd showers..." <<std::endl;
        return false;
       }   

      auto const& shr1 = ev_s->at(0) ;
      auto const& shr2 = ev_s->at(1);

      geoalgo::Vector_t rev_shr1(-1.*shr1.Direction()) ;
      geoalgo::Vector_t rev_shr2(-1.*shr2.Direction()) ;

      // Make the backwards projection for the showers
      auto shr1_bkwrd_hl = ::geoalgo::HalfLine_t(shr1.ShowerStart(),rev_shr1);
      auto shr2_bkwrd_hl = ::geoalgo::HalfLine_t(shr2.ShowerStart(),rev_shr2);

      auto IP = pow(_geoAlgo.SqDist(shr1_bkwrd_hl,shr2_bkwrd_hl),0.5);

      // CCNC the Opening angle of the showers
      double oangle = acos( shr1.Direction().Dot(shr2.Direction())) ;

      // CCNC the vertex point of the two showers. the true designated backwards project
      geoalgo::Point_t vertex(3);

      auto st1 = shr1.ShowerStart();
      auto st2= shr2.ShowerStart();
      auto dir1 = shr1.Direction();
      auto dir2 = shr2.Direction();
      geoalgo::HalfLine_t shr1_hl(st1.X(),st1.Y(),st1.Z(),dir1.X(),dir1.Y(), dir1.Z() );
      geoalgo::HalfLine_t shr2_hl(st2.X(),st2.Y(),st2.Z(),dir2.X(),dir2.Y(), dir2.Z() );

      _geoAlgo.commonOrigin(shr1_hl, shr2_hl, vertex, true);

      // CCNC Diretion of two correlated shower
      geoalgo::Vector_t momentum(3);// need to fill out
      geoalgo::Vector_t mom_vect(shr2.Direction()*shr1.Energy(2) +shr1.Direction()*shr2.Energy(2)) ;
      auto tot_pi0_mom = sqrt(pow(mom_vect[0],2) + pow(mom_vect[1],2) + pow(mom_vect[2],2) );

      //===========================================
      auto radL_shr1 = vertex.Dist(shr1.ShowerStart());
      auto radL_shr2 = vertex.Dist(shr2.ShowerStart());
      //===========================================

      _pi0_mass      = sqrt(2 * shr1.Energy() * shr2.Energy() *(1.-cos(oangle)));
      _pi0_mom       = tot_pi0_mom;
      _pi0_oangle    = oangle;
      _pi0_IP        = IP ;
      _pi0_low_shrE  = shr1.Energy() < shr2.Energy() ? shr1.Energy() : shr2.Energy() ;
      _pi0_high_shrE = shr1.Energy() < shr2.Energy() ? shr2.Energy() : shr1.Energy() ;
      _pi0_low_radL  = shr1.Energy() < shr2.Energy() ? radL_shr1 : radL_shr2 ;
      _pi0_high_radL = shr1.Energy() < shr2.Energy() ? radL_shr2 : radL_shr1 ;

    }
 
    if ( _get_single_shower_info ){
      auto ev_s = storage->get_data<event_shower>("pi0_1gamma_candidate_showers");
      if ( !ev_s ) return false;

      geoalgo::Point_t vertex(3);
      vertex[0] = vtx.X();
      vertex[1] = vtx.Y();
      vertex[2] = vtx.Z();

      _gamma_E = ev_s->at(0).Energy(2); 
      _gamma_RL = vertex.Dist(ev_s->at(0).ShowerStart());

    }

    //std::cout<<"Event and ID: "<<_event<<", "<<_bkgd_v[_bkgd_id]<<"\n\n";
    
    _tree->Fill();    

    return true;
  }

  bool BackgroundBT::finalize() {

    std::cout<<"Signals: "<<std::endl ;
    std::cout<<"Total CCpi0 : "<<_n_cc1pi0<<std::endl; 

    // Note that cc other includes secondary pi0s.
    std::cout<<"\nBackgrounds: "<<std::endl;
    std::cout<<"0) Noise : "<<_n_noise<< std::endl;
    std::cout<<"1) Cosmic : "<<_n_cosmic<< std::endl;
    std::cout<<"2) CC 1pi0 : "<<_n_cc1pi0<<std::endl;
    std::cout<<"3) NC 1pi0 : "<<_n_nc1pi0<<std::endl;
    std::cout<<"4) CC 1pi0 out FV: "<<_n_cc1pi0_outFV<<std::endl;
    std::cout<<"5) CC mult pi0: "<<_n_multpi0<<std::endl;
    std::cout<<"6) Nue: "<<_n_nue<<std::endl;
    std::cout<<"7) Antinumu: "<<_n_antimu<<std::endl;
    std::cout<<"8) CC Cex: "<<_n_cccex<<std::endl;
    std::cout<<"9) NC Cex: "<<_n_nccex<<std::endl;
    std::cout<<"10) CC Gamma : "<<_n_gamma<<std::endl;
    std::cout<<"11) Kaon Decay: "<<_n_kaondecay<<std::endl;
    std::cout<<"12) CC Other : "<<_n_ccother<<std::endl;
    std::cout<<"13) NC Other : "<<_n_ncother<<std::endl;
    std::cout<<"14) Other: "<<_n_other<< std::endl;

    //std::cout<<"Total accounted backgrounds: "<< _n_other + _n_cosmic + _n_nc1pi0 + _n_nc0pi0 + _n_cc0pi0 <<std::endl ;

    for (int i = 0; i < _event_list.size(); i++){
      std::cout<<_event_list[i]<<std::endl;
      }

    if ( _fout ){
      _fout->cd();
      _tree->Write();
    }
  
    return true;
  }

}
#endif