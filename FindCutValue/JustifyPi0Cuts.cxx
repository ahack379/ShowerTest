#ifndef LARLITE_JUSTIFYCUTS_CXX
#define LARLITE_JUSTIFYCUTS_CXX

#include "JustifyPi0Cuts.h"
#include "DataFormat/shower.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/cluster.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"

#include "DataFormat/hit.h"

namespace larlite {

  bool JustifyPi0Cuts::initialize() {

    _n_other = 0;    // 0 
    _n_cosmic = 0;   // 1
    _n_cc1pi0 = 0;   // 2 
    _n_cc0pi0 = 0;   // 3
    _n_nc1pi0 = 0;   // 4 
    _n_nc0pi0 = 0;   // 5

    if( !_gamma_tree ){
      _gamma_tree = new TTree("gamma_tree","");
      _gamma_tree->Branch("_event",&_event,"event/I");
      _gamma_tree->Branch("bkgd_id",&_bkgd_id,"bkgd_id/I");
      _gamma_tree->Branch("_gamma_low_E",&_gamma_low_E,"gamma_low_E/F");
      _gamma_tree->Branch("_gamma_high_E",&_gamma_high_E,"gamma_high_E/F");
      _gamma_tree->Branch("_gamma_low_RL",&_gamma_low_RL,"gamma_low_RL/F");
      _gamma_tree->Branch("_gamma_high_RL",&_gamma_high_RL,"gamma_high_RL/F");
      _gamma_tree->Branch("_gamma_oangle",&_gamma_oangle,"gamma_oangle/F");
      _gamma_tree->Branch("_gamma_IP",&_gamma_IP,"gamma_IP/F");
      _gamma_tree->Branch("_gamma_low_matched",&_gamma_low_matched,"gamma_low_matched/B");
      _gamma_tree->Branch("_gamma_high_matched",&_gamma_high_matched,"gamma_high_matched/B");
      _gamma_tree->Branch("_gamma1_vtx_IP",&_gamma1_vtx_IP,"gamma1_vtx_IP/F");
      _gamma_tree->Branch("_gamma2_vtx_IP",&_gamma2_vtx_IP,"gamma2_vtx_IP/F");
      _gamma_tree->Branch("_pi0_mass",&_pi0_mass,"pi0_mass/F");
      _gamma_tree->Branch("_pi0_mom",&_pi0_mom,"pi0_mom/F");
      _gamma_tree->Branch("_event_type",&_event_type,"event_type/I");
      _gamma_tree->Branch("_nu_pdg",&_nu_pdg,"nu_pdg/I");
      _gamma_tree->Branch("_isCC",&_isCC,"isCC/B");
      _gamma_tree->Branch("_found_pi0",&_found_pi0,"found_pi0/B");

      _gamma_tree->Branch("_pi0_origin",&_pi0_origin,"pi0_origin/B");
      _gamma_tree->Branch("_pi0_type",&_pi0_type,"pi0_type/B");

      _gamma_tree->Branch("_n_nu_origin_pi0",&_n_nu_origin_pi0,"n_nu_origin_pi0/I");
      _gamma_tree->Branch("_gamma_low_purity",&_gamma_low_purity,"gamma_low_purity/F");
      _gamma_tree->Branch("_gamma_low_complete",&_gamma_low_complete,"gamma_low_complete/F");
      _gamma_tree->Branch("_gamma_high_complete",&_gamma_high_complete,"gamma_high_complete/F");
      _gamma_tree->Branch("_gamma_high_purity",&_gamma_high_purity,"gamma_high_purity/F");

      _gamma_tree->Branch("_mc_startx",&_mc_startx,"mc_startx/F");
      _gamma_tree->Branch("_mc_starty",&_mc_starty,"mc_starty/F");
      _gamma_tree->Branch("_mc_startz",&_mc_startz,"mc_startz/F");
      _gamma_tree->Branch("_gamma_startx",&_gamma_startx,"gamma_startx/F");
      _gamma_tree->Branch("_gamma_starty",&_gamma_starty,"gamma_starty/F");
      _gamma_tree->Branch("_gamma_startz",&_gamma_startz,"gamma_startz/F");

      _gamma_tree->Branch("_mc_low_dirx",&_mc_low_dirx,"mc_low_dirx/F");
      _gamma_tree->Branch("_mc_low_diry",&_mc_low_diry,"mc_low_diry/F");
      _gamma_tree->Branch("_mc_low_dirz",&_mc_low_dirz,"mc_low_dirz/F");
      _gamma_tree->Branch("_mc_high_dirx",&_mc_high_dirx,"mc_high_dirx/F");
      _gamma_tree->Branch("_mc_high_diry",&_mc_high_diry,"mc_high_diry/F");
      _gamma_tree->Branch("_mc_high_dirz",&_mc_high_dirz,"mc_high_dirz/F");
      _gamma_tree->Branch("_gamma_low_dirx",&_gamma_low_dirx,"gamma_low_dirx/F");
      _gamma_tree->Branch("_gamma_low_diry",&_gamma_low_diry,"gamma_low_diry/F");
      _gamma_tree->Branch("_gamma_low_dirz",&_gamma_low_dirz,"gamma_low_dirz/F");
      _gamma_tree->Branch("_gamma_high_dirx",&_gamma_high_dirx,"gamma_high_dirx/F");
      _gamma_tree->Branch("_gamma_high_diry",&_gamma_high_diry,"gamma_high_diry/F");
      _gamma_tree->Branch("_gamma_high_dirz",&_gamma_high_dirz,"gamma_high_dirz/F");
      _gamma_tree->Branch("_res_high",&_res_high,"res_high/F");
      _gamma_tree->Branch("_res_low",&_res_low,"res_low/F");

      _gamma_tree->Branch("_mc_low_startx",&_mc_low_startx,"mc_low_startx/F");
      _gamma_tree->Branch("_mc_low_starty",&_mc_low_starty,"mc_low_starty/F");
      _gamma_tree->Branch("_mc_low_startz",&_mc_low_startz,"mc_low_startz/F");
      _gamma_tree->Branch("_gamma_low_startx",&_gamma_low_startx,"gamma_low_startx/F");
      _gamma_tree->Branch("_gamma_low_starty",&_gamma_low_starty,"gamma_low_starty/F");
      _gamma_tree->Branch("_gamma_low_startz",&_gamma_low_startz,"gamma_low_startz/F");
      _gamma_tree->Branch("_mc_high_startx",&_mc_high_startx,"mc_high_startx/F");
      _gamma_tree->Branch("_mc_high_starty",&_mc_high_starty,"mc_high_starty/F");
      _gamma_tree->Branch("_mc_high_startz",&_mc_high_startz,"mc_high_startz/F");
      _gamma_tree->Branch("_gamma_high_startx",&_gamma_high_startx,"gamma_high_startx/F");
      _gamma_tree->Branch("_gamma_high_starty",&_gamma_high_starty,"gamma_high_starty/F");
      _gamma_tree->Branch("_gamma_high_startz",&_gamma_high_startz,"gamma_high_startz/F");


    }

    if( !_one_gamma_tree ){
      _one_gamma_tree = new TTree("one_gamma_tree","");
      _one_gamma_tree->Branch("bkgd_id",&_bkgd_id,"bkgd_id/I");
      _one_gamma_tree->Branch("_event",&_event,"event/I");
      _one_gamma_tree->Branch("_gamma_E",&_gamma_E,"gamma_E/F");
      _one_gamma_tree->Branch("_gamma_RL",&_gamma_RL,"gamma_RL/F");
      _one_gamma_tree->Branch("_gamma_vtx_IP",&_gamma_vtx_IP,"gamma_vtx_IP/F");
      _one_gamma_tree->Branch("_gamma_matched",&_gamma_matched,"gamma_matched/B");
      _one_gamma_tree->Branch("_event_type",&_event_type,"event_type/I");
      _one_gamma_tree->Branch("_nu_pdg",&_nu_pdg,"nu_pdg/I");
      _one_gamma_tree->Branch("_isCC",&_isCC,"isCC/B");
      _one_gamma_tree->Branch("_found_pi0",&_found_pi0,"found_pi0/B");
      _one_gamma_tree->Branch("_n_nu_origin_pi0",&_n_nu_origin_pi0,"n_nu_origin_pi0/I");
      _one_gamma_tree->Branch("_gamma_purity",&_gamma_purity,"gamma_purity/F");
      _one_gamma_tree->Branch("_gamma_complete",&_gamma_complete,"gamma_complete/F");

      _one_gamma_tree->Branch("_mc_dirx",&_mc_dirx,"mc_dirx/F");
      _one_gamma_tree->Branch("_mc_diry",&_mc_diry,"mc_diry/F");
      _one_gamma_tree->Branch("_mc_dirz",&_mc_dirz,"mc_dirz/F");
      _one_gamma_tree->Branch("_gamma_dirx",&_gamma_dirx,"gamma_dirx/F");
      _one_gamma_tree->Branch("_gamma_diry",&_gamma_diry,"gamma_diry/F");
      _one_gamma_tree->Branch("_gamma_dirz",&_gamma_dirz,"gamma_dirz/F");

      _one_gamma_tree->Branch("_res",&_res,"res/F");
      _one_gamma_tree->Branch("_gamma_origin",&_gamma_origin,"gamma_origin/I");
      _one_gamma_tree->Branch("_gamma_type",&_gamma_type,"gamma_type/I");
    }

    if(!_tree){
      _tree = new TTree("tree","");
      _tree->Branch("_mu_startx",&_mu_startx,"mu_startx/F");
      _tree->Branch("_mu_starty",&_mu_starty,"mu_starty/F");
      _tree->Branch("_mu_startz",&_mu_startz,"mu_startz/F");
      _tree->Branch("_mu_endx",&_mu_endx,"mu_endx/F");
      _tree->Branch("_mu_endy",&_mu_endy,"mu_endy/F");
      _tree->Branch("_mu_endz",&_mu_endz,"mu_endz/F");
      _tree->Branch("_mu_mom",&_mu_mom,"mu_mom/F");
      _tree->Branch("_mu_len",&_mu_len,"mu_len/F");
      _tree->Branch("_mu_angle",&_mu_angle,"mu_angle/F");
      _tree->Branch("_mu_phi",&_mu_phi,"mu_phi/F");
      _tree->Branch("_mult",&_mult,"mult/F");
  
      _tree->Branch("bkgd_id",&_bkgd_id,"bkgd_id/I");
      _tree->Branch("nshrs",&_nshrs,"nshrs/I");
    }

    if(!_compare_tree){
      _compare_tree = new TTree("compare_tree","");
      _compare_tree->Branch("_event",&_event,"event/I");
      _compare_tree->Branch("bkgd_id",&_bkgd_id,"bkgd_id/I");
      _compare_tree->Branch("_reco_E",&_reco_E,"reco_E/F");
      _compare_tree->Branch("_reco_startx",&_reco_startx,"reco_startx/F");
      _compare_tree->Branch("_reco_starty",&_reco_starty,"reco_starty/F");
      _compare_tree->Branch("_reco_startz",&_reco_startz,"reco_startz/F");
      _compare_tree->Branch("_reco_start3D",&_reco_start3D,"reco_start3D/F");
      _compare_tree->Branch("_reco_dirx",&_reco_dirx,"reco_dirx/F");
      _compare_tree->Branch("_reco_diry",&_reco_diry,"reco_diry/F");
      _compare_tree->Branch("_reco_dirz",&_reco_dirz,"reco_dirz/F");
      _compare_tree->Branch("_reco_dot",&_reco_dot,"reco_dot/F");

      _compare_tree->Branch("_mc_E",&_mc_E,"mc_E/F");
      _compare_tree->Branch("_mc_startx",&_mc_startx,"mc_startx/F");
      _compare_tree->Branch("_mc_starty",&_mc_starty,"mc_starty/F");
      _compare_tree->Branch("_mc_startz",&_mc_startz,"mc_startz/F");
      _compare_tree->Branch("_mc_dirx",&_mc_dirx,"mc_dirx/F");
      _compare_tree->Branch("_mc_diry",&_mc_diry,"mc_diry/F");
      _compare_tree->Branch("_mc_dirz",&_mc_dirz,"mc_dirz/F");
    }

    _SCE = new larutil::SpaceChargeMicroBooNE();
    _time2cm = larutil::GeometryHelper::GetME()->TimeToCm();

    _event = -1;

    return true;
  }

  void JustifyPi0Cuts::clear(){
  
    _pi0_mass      = -10;
    _pi0_mom       = -10;
    _gamma_oangle    = -10;
    _gamma_low_E  = -10;
    _gamma_high_E = -10;
    _gamma_low_RL  = -10;
    _gamma_high_RL = -10;
    _gamma_IP = -10;
    _gamma1_vtx_IP = -10;
    _gamma2_vtx_IP = -10;
    _gamma_low_matched = false;
    _gamma_high_matched = false;

    _gamma_E = -10;
    _gamma_RL = -10;
    _gamma_vtx_IP = -10;
    _gamma_matched =false;

    _gamma_low_purity = -999;
    _gamma_low_complete = -999;
    _gamma_high_purity = -999;
    _gamma_high_complete = -999;

    _gamma_low_dirx = -999 ;
    _gamma_low_diry = -999 ;
    _gamma_low_dirz = -999 ;
    _gamma_high_dirx = -999 ;
    _gamma_high_diry = -999 ;
    _gamma_high_dirz = -999 ;
    _mc_low_dirx = -999 ;
    _mc_low_diry = -999 ;
    _mc_low_dirz = -999 ;
    _mc_high_dirx = -999 ;
    _mc_high_diry = -999 ;
    _mc_high_dirz = -999 ;

    _mc_startx = -999;
    _mc_starty = -999;
    _mc_startz = -999;
    _gamma_startx = -999;
    _gamma_starty = -999;
    _gamma_startz = -999;

    _mc_low_startx = -999;
    _mc_low_starty = -999;
    _mc_low_startz = -999;
    _mc_high_startx = -999;
    _mc_high_starty = -999;
    _mc_high_startz = -999;
    _gamma_low_startx = -999;
    _gamma_low_starty = -999;
    _gamma_low_startz = -999;
    _gamma_high_startx = -999;
    _gamma_high_starty = -999;
    _gamma_high_startz = -999;

    _res_low = -999;
    _res_high = -999;
    _res = -999;

    _pi0_origin = false;
    _pi0_type = false;

    _event_type = -1;

    _event_type = -1;

    _gamma_purity = -999;
    _gamma_complete = -999;
    _gamma_origin = -999;
    _gamma_type = -999;

    _mc_dirx = -999;
    _mc_diry = -999;
    _mc_dirz = -999;
    _gamma_dirx = -999 ;
    _gamma_diry = -999 ;
    _gamma_dirz = -999 ;

    _event_type = -1;

    _mu_startx     = -1000;
    _mu_starty     = -1000;
    _mu_startz     = -1000;
    _mu_endx     = -1000;
    _mu_endy     = -1000;
    _mu_endz     = -1000;
    _mu_len        = -10;
    _mu_len        = -10;
    _mu_mom        = -10;
    _mu_angle      = -10;
    _mu_phi        = -10;
    _mult          = 0;
    
    _bkgd_id       = -1 ;
    _nshrs         = -1;
  }

  void JustifyPi0Cuts::clearCompare(){
  
    _reco_E = -999;
    _reco_dot = -999;
    _reco_start3D = -999;
    _reco_startx = -999;
    _reco_starty = -999;
    _reco_startz = -999;
    _reco_dirx = -999;
    _reco_diry = -999;
    _reco_dirz = -999;

    _mc_startx = -999;
    _mc_starty = -999;
    _mc_startz = -999;
    _mc_dirx = -999;
    _mc_diry = -999;
    _mc_dirz = -999;
  
  }
  
  bool JustifyPi0Cuts::analyze(storage_manager* storage) {

    _event++;
    clear();

    auto ev_s = storage->get_data<event_shower>("showerreco");
    if( !ev_s || !ev_s->size() || ev_s->size() < 1 )
      return false;

    std::cout<<"\n\nEvent : "<<_event <<std::endl;

      auto ev_t = storage->get_data<event_track>("numuCC_track");

      if( !ev_t || !ev_t->size() ){
        std::cout<<"No tagged track; what??" <<std::endl;
        return false;
      }

      auto ev_v = storage->get_data<event_vertex>("numuCC_vertex");

      if( !ev_v || !ev_v->size() ){
        std::cout<<"No tagged vertex; what??" <<std::endl;
        return false;
      }

      auto t = ev_t->at(0) ;
      auto v = ev_v->at(0) ;

      _mu_startx = t.Vertex().X() ;
      _mu_starty = t.Vertex().Y() ;
      _mu_startz = t.Vertex().Z() ;
      _mu_endx = t.End().X() ;
      _mu_endy = t.End().Y() ;
      _mu_endz = t.End().Z() ;
      //_mu_mom  = t.VertexMomentum() ;
      _mu_len  = t.Length(0); // Calculates the length from point 0 to end
      _mu_angle = t.VertexDirection().Theta() ;
      _mu_phi = t.VertexDirection().Phi() ; 

      std::vector<double> dir = { (_mu_endx - _mu_startx) / _mu_len,
                                  (_mu_endy - _mu_starty) / _mu_len,
                                  (_mu_endz - _mu_startz) / _mu_len };

      auto dir_start = t.VertexDirection();
      std::vector<double> other_dir = { dir_start.X(), dir_start.Y(), dir_start.Z() };  

      float dotProd = dir.at(0) * other_dir.at(0) + dir.at(1) * other_dir.at(1) +  dir.at(2) * other_dir.at(2) ;

      if( dotProd < 0 ) { 
         TVector3 new_dir(-dir_start.X(),-dir_start.Y(),-dir_start.Z());
         _mu_angle = new_dir.Theta();
         _mu_phi = new_dir.Phi();
      }   

      auto ev_trk = storage->get_data<event_track>("pandoraNu");

      for ( auto const & t : *ev_trk ){
        auto st = t.Vertex() ;
        auto end = t.Vertex() ;

        auto dist_st = sqrt( pow(st.X() - v.X(),2) + pow(st.Y() - v.Y(),2) + pow(st.Z() - v.Z(),2) );
        auto dist_end = sqrt( pow(end.X() - v.X(),2) + pow(end.Y() - v.Y(),2) + pow(end.Z() - v.Z(),2) );

        if (dist_st < 3 || dist_end < 3)
          _mult ++ ;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////
      // Want to be able to access the origin of the tagged muon. Thus, need to find it, and 
      // Ask for its origin.  Need to match to MCtrack to do this
      auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
      if(!ev_mctruth || !ev_mctruth->size() ) {
        std::cout<<"Event has no mctruth info "<<std::endl;
        return false;
      }

      auto ev_mctrk = storage->get_data<event_mctrack>("mcreco");
      if ( !ev_mctrk || !ev_mctrk->size() ) {std::cout<<"No MCTrack!" <<std::endl ; return false; }


      auto & truth = ev_mctruth->at(0);
      auto & nu  = truth.GetNeutrino();

      double xyz[3] = {0.};
      auto traj = nu.Nu().Trajectory();
      xyz[0] = traj.at(traj.size() - 1).X();
      xyz[1] = traj.at(traj.size() - 1).Y();
      xyz[2] = traj.at(traj.size() - 1).Z();

      // Correct for space charge effect
      auto tvtx = traj.at(traj.size() - 1).T(); // ns
      auto vtxtick = (tvtx / 1000.) * 2.; // time in tick :
      //auto vtxtimecm = vtxtick * _time2cm; // time in cm :

      auto ev_mcc = storage->get_data<event_cluster>("mccluster");

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

      auto ev_hit = storage->get_data<larlite::event_hit>("gaushit");

      if (!ev_hit || ev_hit->size() == 0){
        std::cout << "No hits! exit" << std::endl;
        return false;
      }      

    // Fill multiplicity info + find ID of pandora track that is numuCC_track
    int min_trk_dist = 1e9;
    int min_trk_dist_it = -1;

    for ( int ii = 0; ii < ev_trk->size(); ii++){

        auto t = ev_trk->at(ii);
        auto st = t.Vertex() ;
        auto end = t.End() ;

        auto tag_end = ev_t->at(0).End() ;
        auto dist = sqrt( pow(tag_end.X() - end.X(),2) + pow(tag_end.Y() - end.Y(),2) + pow(tag_end.Z() - end.Z(),2) );
        if ( dist < min_trk_dist ){
          min_trk_dist = dist ;
          min_trk_dist_it = ii ;
        }
    }

      std::vector<int> pur_ctr_v ;
      std::vector<float> cw_pur_ctr_v ;

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

      // for each reco cluster, find the origin of all hits and calc purity/completeness 
      // the "...size()+1" is to account for noise category
      pur_ctr_v.resize(ass_mcclus_v.size()+1,0) ;
      cw_pur_ctr_v.resize(ass_mcclus_v.size()+1,0) ;

      int max_hits = -1;
      int max_cw_hits = -1;
      int max_cid = -1 ;
      float tot_reco_cw_hits = 0;

      // This particular headache is necessary because there are no gaushit associations
      // to pandoraNu tracks, only pandoraNuCosmicRemoval associations. 
      // We're using gaushit for everything else including clustering + mccluster building though,
      // so need to find the gaushits that correspond to the cosmic removal hits 
      // and store these IDs for backtracker to work.
      std::vector<int> tag_trk_gaushit_v;
      for(int i = 0; i < ass_hit_v.at(min_trk_dist_it).size(); i++){
         auto hid = ass_hit_v.at(min_trk_dist_it).at(i) ;
         auto h = ev_hit_cosRem->at(hid);
         if ( h.WireID().Plane != 2 ) continue;

         for(int j = 0; j < ev_hit->size(); j++){
           auto hj = ev_hit->at(j) ;
           if ( hj.WireID().Plane != 2 ) continue;

           if ( hj.PeakTime() == h.PeakTime() && hj.WireID().Wire == h.WireID().Wire )
             tag_trk_gaushit_v.emplace_back(j);
         }
      }

      //std::cout<<"Tagged track : "<<tag_trk_gaushit_v.size()<<std::endl ;

      // Now calculate purity and completeness for muon track
      for(int i = 0; i < tag_trk_gaushit_v.size(); i++){
         auto hid = tag_trk_gaushit_v.at(i) ;
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
         // our "else" here is to account for noise category in which there is no corresponding true charge
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

         auto tot_mc_hits =  ass_mcclus_v[max_cid].size();
         auto tot_reco_hits = tag_trk_gaushit_v.size() ; //ass_hit_v.at(min_trk_dist_it).size();         

       }
      else 
         _bkgd_id = 0 ;

      ///////////////////////////////////////////////////////////////////////////////////////////////
      /// Now count the number of backgrounds and signals
      ///////////////////////////////////////////////////////////////////////////////////////////////
      if( _bkgd_id == -1 ){
        auto ev_mcs = storage->get_data<event_mcshower>("mcreco") ;
        if ( !ev_mcs || !ev_mcs->size() ) {std::cout<<"No MCShower!" <<std::endl ; return false; }

        bool infv = true;
        if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1]> 96.5 || xyz[1] < -96.5 || xyz[2]< 10 || xyz[2]> 1026.8 )
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

        _nu_pdg = nu.Nu().PdgCode() ;
        _isCC = nu.CCNC() == 0 ? true : false ;
        _n_nu_origin_pi0 = n_pi0 ;

       auto mcclus = ev_mcc->at(max_cid) ;
       auto mu_origin = mcclus.Width() ;

       //std::cout<<"muon purity : "<<_mu_purity<<", "<<_mu_complete<<", "<<_mu_origin<<", "<<_mu_type<<std::endl ; 

       // Remember that I've repurposed the cluster width variable to store info 
       // about whether this mccluster (which we've identified as the match to 
       // our reco particle at this stage) is neutrino or cosmic in origin 

        //std::cout<<"mccluster width : "<<mcclus.Width()<<std::endl ;
        if( mu_origin == 2)
          _bkgd_id = 1;
        else if( abs(nu.Nu().PdgCode()) == 12 ){
          _bkgd_id = 6 ;
        }
        else if( nu.Nu().PdgCode() == -14 ){
          _bkgd_id = 7 ;
        }
        // frmo here we can assume we have a muon neutrino
        else if( nu.Nu().PdgCode() == 14 && n_pi0 == 1 && nu.CCNC() == 0 && infv){
          _bkgd_id = 2;
        }
        else if( nu.CCNC() == 1 && n_pi0 > 0 ) {
          _bkgd_id = 3;
        }
        else if ( nu.CCNC() == 0 && n_pi0 == 1 && !infv )
          _bkgd_id = 4 ;
       else if( nu.CCNC() == 0 && n_pi0 > 1 ) {
          _bkgd_id = 5;
        }
        else if( n_pi0 == 0 && n_gamma > 0 ){
          _bkgd_id = 10 ;
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
          }
          else if( charge_ex && nu.CCNC() == 1 ) {
            _bkgd_id = 9;
          }
          else if( kaon_decay ){
            _bkgd_id = 11;
          }
          else if( !charge_ex && nu.CCNC() == 0 ){
            _bkgd_id = 12;

          }
          else if( !charge_ex && nu.CCNC() == 1 ){
            _bkgd_id = 13;
          }
          else {
            _bkgd_id = 14;
          }
        }
      }

    //auto ev_s = storage->get_data<event_shower>("showerreco");
    _nshrs = ev_s->size() ;
   
    _tree->Fill();

    //if( !ev_s || !ev_s->size() || ev_s->size() < 1 ){
    //  std::cout<<"Not enough reco'd showers..." <<ev_s->size()<<std::endl;
    //  return false;
    //} 

    auto ev_mcshr = storage->get_data<event_mcshower>("mcreco");

    // Now deal with shower comparisons
    std::vector<int> shr_ids;

    bool found_pi0 = false;
    int pi0_ctr = 0;

    for ( int si = 0; si < ev_mcshr->size(); si++){

      auto s = ev_mcshr->at(si);
      auto st = s.Start();
      auto dist = sqrt( pow(st.X() - xyz[0],2) + pow(st.Y() - xyz[1],2) + pow(st.Z() - xyz[2],2) );

      if ( s.DetProfile().E() > 0 && s.MotherPdgCode() == 111 && s.Origin() == 1 ){ found_pi0 = true; }

      if (found_pi0 && dist < 0.00001 && nu.Nu().PdgCode() == 14 && nu.CCNC() == 0 ){
        pi0_ctr ++ ;
        shr_ids.emplace_back(si) ;
        _event_type = 0;
      }
    }
    
    _found_pi0 = found_pi0 ;

    if ( pi0_ctr > 2 && shr_ids.size() > 0 )  _event_type = 1;  

    if ( found_pi0 && shr_ids.size() == 0 ){
      for ( int si = 0; si < ev_mcshr->size(); si++){

        auto s = ev_mcshr->at(si);
        auto st = s.Start();
        auto dist = sqrt( pow(st.X() - xyz[0],2) + pow(st.Y() - xyz[1],2) + pow(st.Z() - xyz[2],2) );

        if ( s.DetProfile().E() > 0 && s.MotherPdgCode() == 111 ){
          shr_ids.emplace_back(si) ;
	      _event_type = 1;
	    }
	  }
    }
      
    if ( _event_type == -1 ) _event_type = 2;

    std::multimap<float,std::pair<int,int>> mc_reco_map ;    

    //std::cout<<" shr ids: "<<shr_ids.size()<<std::endl;
    // Match showers
    for( auto const & mc_id : shr_ids ){
      auto mcs_i = ev_mcshr->at(mc_id);
      auto mag_mcs = sqrt( pow(mcs_i.DetProfile().Px(),2) + pow(mcs_i.DetProfile().Py(),2) + pow(mcs_i.DetProfile().Pz(),2) );

      for( int reco_id = 0; reco_id < ev_s->size(); reco_id++ ){

        auto recos_i = ev_s->at(reco_id) ;

        auto mag_reco = sqrt( pow(recos_i.Direction().Px(),2) + pow(recos_i.Direction().Py(),2) + pow(recos_i.Direction().Pz(),2) );
        auto dot = mcs_i.DetProfile().Px() * recos_i.Direction().Px() +
                   mcs_i.DetProfile().Py() * recos_i.Direction().Py() +
                   mcs_i.DetProfile().Pz() * recos_i.Direction().Pz() ;
        dot /= ( mag_mcs * mag_reco );
        //std::cout<<" dots and things: "<<dot<<std::endl ;

        //if ( fabs(dot) > 1 ) std::cout<<"DOT ! " <<dot <<std::endl ;

        mc_reco_map.emplace(fabs(1./dot),std::make_pair(mc_id,reco_id)) ;
        //std::cout<<"mc and reco ids: "<<mc_id<<", "<<reco_id <<", "<<ev_s->size()<<std::endl ;
      }
    } 

    int reco_g1_id = -1, reco_g2_id = -1;
    int mc_g1_id = -1, mc_g2_id = -1;
    float dot_g1 = -10, dot_g2 = -10;

    for ( auto const & m : mc_reco_map ){

      auto score = 1./m.first ;
       
      //std::cout<<"Some scores are..." <<score<<std::endl;
    
      //if ( score < 0.9 ) break; 

      if ( reco_g1_id == -1 ){
        mc_g1_id   = m.second.first ;
        reco_g1_id = m.second.second ;
        dot_g1 = score ;
        if ( ev_s->size() == 1 ) break;
      }
      else{
        //if ( m.second.first == mc_g1_id || m.second.second == reco_g1_id ) continue;
        
        mc_g2_id = m.second.first;       
        reco_g2_id = m.second.second ;       
        dot_g2 = score ;
        break;
      }
    }
   // if ( mc_g2_id != -1 ) {

   // std::cout<<"Mc and reco ids : "<<mc_g1_id<<", "<<reco_g1_id<<", "<<dot_g1<<", energy: "<<ev_s->at(reco_g1_id).Energy(2)<<", "<<ev_mcshr->at(mc_g1_id).DetProfile().E()<<std::endl;
   // std::cout<<"Mc and reco ids : "<<mc_g2_id<<", "<<reco_g2_id<<", "<<dot_g2<<", energy: "<<ev_s->at(reco_g2_id).Energy(2)<<", "<<ev_mcshr->at(mc_g2_id).DetProfile().E()<<std::endl;

   //}

    // Get the association from cluster -> hit
    auto const & ev_clus = storage->get_data<event_cluster>("ImageClusterHit");
    auto const & ev_ass_c = storage->get_data<larlite::event_ass>("ImageClusterHit");
    auto const & ass_imageclus_v = ev_ass_c->association(ev_clus->id(), ev_hit->id());

    //auto ev_mcs = storage->get_data<event_mcshower>("mcreco") ;
    //if ( !ev_mcs || !ev_mcs->size() ) {std::cout<<"No MCShower!" <<std::endl ; return false; }
    //ev_mcshr

    // Get the association from shower -> cluster
    auto ev_ass_s = storage->get_data<larlite::event_ass>("showerreco");
    auto const& ass_showerreco_v = ev_ass_s->association(ev_s->id(), ev_clus->id());
    std::cout<<"Background : "<<_event<<", "<<_bkgd_id <<std::endl ;

    for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

       clearCompare();

       auto const& shr1 = ev_s->at(s1);

       _gamma_E  = shr1.Energy(2)  ;
       ::geoalgo::Point_t vertex_reco(v.X(),v.Y(),v.Z());
        _gamma_RL = vertex_reco.Dist(shr1.ShowerStart());

       geoalgo::Vector_t rev_shr(-1.*shr1.Direction()) ;
       auto shr_bkwrd_hl = ::geoalgo::HalfLine_t(shr1.ShowerStart(),rev_shr);
       _gamma_vtx_IP = _geoAlgo.SqDist(vertex_reco, shr_bkwrd_hl) ;

       auto dir1 = shr1.Direction();
       _gamma_dirx = dir1.X();
       _gamma_diry = dir1.Y();
       _gamma_dirz = dir1.Z();

       if ( s1 == reco_g1_id ){

	     auto mcs = ev_mcshr->at(mc_g1_id);
	     _res = dot_g1 ;
         _mc_dirx = mcs.DetProfile().Px();  
         _mc_diry = mcs.DetProfile().Py();  
         _mc_dirz = mcs.DetProfile().Pz();  
         //_mc_startx = mcs.DetProfile().X();
         //_mc_starty = mcs.DetProfile().Y();
         //_mc_startz = mcs.DetProfile().Z();

         _gamma_startx = shr1.ShowerStart().X(); 
         _gamma_starty = shr1.ShowerStart().Y(); 
         _gamma_startz = shr1.ShowerStart().Z(); 

         auto mcx = mcs.DetProfile().X() ;
         auto mcy = mcs.DetProfile().Y() ;
         auto mcz = mcs.DetProfile().Z() ;
         auto mct = mcs.DetProfile().T() ;
         auto vtxtick = (mct/ 1000.) * 2.; 
         auto vtxtimecm = vtxtick * _time2cm; 
         auto sce_corr = _SCE->GetPosOffsets(mcx,mcy,mcz);

         _mc_startx = mcx + vtxtimecm + 0.7 - sce_corr.at(0);
         _mc_starty = mcy + sce_corr.at(1);
         _mc_startz = mcz + sce_corr.at(2);
       }

       if ( s1 == reco_g2_id ){
	     auto mcs = ev_mcshr->at(mc_g2_id);
	     _res = dot_g2 ;
         _mc_dirx = mcs.DetProfile().Px();  
         _mc_diry = mcs.DetProfile().Py();  
         _mc_dirz = mcs.DetProfile().Pz();  
         _gamma_startx = shr1.ShowerStart().X(); 
         _gamma_starty = shr1.ShowerStart().Y(); 
         _gamma_startz = shr1.ShowerStart().Z(); 

         auto mcx = mcs.DetProfile().X() ;
         auto mcy = mcs.DetProfile().Y() ;
         auto mcz = mcs.DetProfile().Z() ;
         auto mct = mcs.DetProfile().T() ;
         auto vtxtick = (mct/ 1000.) * 2.; 
         auto vtxtimecm = vtxtick * _time2cm; 
         auto sce_corr = _SCE->GetPosOffsets(mcx,mcy,mcz);

         _mc_startx = mcx + vtxtimecm + 0.7 - sce_corr.at(0);
         _mc_starty = mcy + sce_corr.at(1);
         _mc_startz = mcz + sce_corr.at(2);

       }


       float mc_clus_e = 0.;
       int closest_mcs_id = -1 ;
       float closest_e = 1e9 ;

       int shr1_origin = -1;
       int shr1_type = -1;

       for (size_t j = 0; j < ass_showerreco_v.at(s1).size(); j++ ){

         auto clus_id = ass_showerreco_v.at(s1).at(j);
         auto iclus = ev_clus->at(clus_id);

         int plane = iclus.Plane().Plane ;
         if ( plane != 2 ) continue;

         pur_ctr_v.clear();
         pur_ctr_v.resize(ass_mcclus_v.size()+1,0) ;

         max_hits = -1;
         max_cid = -1 ;

         // Loop through all hits associared to the cluster 
         for ( int k = 0; k < ass_imageclus_v.at(clus_id).size(); k++ ){

           auto hid = ass_imageclus_v.at(clus_id).at(k) ;
           auto h = ev_hit->at(hid);

           if ( _mc_hit_map.find(hid) != _mc_hit_map.end() ){

             auto mcclus_id = _mc_hit_map[hid] ;
             pur_ctr_v[mcclus_id]++ ;

             if( pur_ctr_v[ mcclus_id] > max_hits ){
               max_hits = pur_ctr_v[mcclus_id];
               max_cid = mcclus_id ;
             }
           }   
           else {
             auto mcclus_id = ass_mcclus_v.size() ;
             pur_ctr_v[mcclus_id]++ ; 
             if( pur_ctr_v[mcclus_id] > max_hits ){
               max_hits = pur_ctr_v[mcclus_id];
               max_cid = mcclus_id ; 
             }
           }
         }


          if ( max_cid != ass_mcclus_v.size() && max_cid != -1 ){

            auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
            auto tot_reco_hits = ass_imageclus_v[clus_id].size();

            auto mcclus1 = ev_mcc->at(max_cid) ;
            shr1_origin = mcclus1.Width() ;
            shr1_type = mcclus1.StartOpeningAngle() ;

			_gamma_origin = shr1_origin;
            _gamma_type = shr1_type;
    
            _gamma_purity   = float(max_hits) / tot_reco_hits ;
            _gamma_complete = float(max_hits) / tot_mc_hits ;
          }
       }

       _one_gamma_tree->Fill();

       for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

            if (s2 <= s1) continue;

            //if ( s2 == reco_g2_id ) _gamma_low_matched = true ;

            auto const& shr2 = ev_s->at(s2);

            geoalgo::Vector_t rev_shr1(-1.*shr1.Direction()) ;
            geoalgo::Vector_t rev_shr2(-1.*shr2.Direction()) ;

            // Make the backwards projection for the showers
            auto shr1_bkwrd_hl = ::geoalgo::HalfLine_t(shr1.ShowerStart(),rev_shr1);
            auto shr2_bkwrd_hl = ::geoalgo::HalfLine_t(shr2.ShowerStart(),rev_shr2);

            // Calc the Opening angle of the showers
            double oangle = acos( shr1.Direction().Dot(shr2.Direction())) ;

            // Calc the vertex point of the two showers. the true designated backwards project
            geoalgo::Point_t vertex(3);

            auto st1 = shr1.ShowerStart();
            auto st2 = shr2.ShowerStart();
            auto dir1 = shr1.Direction();
            auto dir2 = shr2.Direction();
            geoalgo::HalfLine_t shr1_hl(st1.X(),st1.Y(),st1.Z(),dir1.X(),dir1.Y(), dir1.Z() );
            geoalgo::HalfLine_t shr2_hl(st2.X(),st2.Y(),st2.Z(),dir2.X(),dir2.Y(), dir2.Z() );

            _geoAlgo.commonOrigin(shr1_hl, shr2_hl, vertex, true);

            // Calc Diretion of two correlated shower
            geoalgo::Vector_t momentum(3);// need to fill out
            geoalgo::Vector_t mom_vect(shr2.Direction()*shr1.Energy(2) +shr1.Direction()*shr2.Energy(2)) ;

            auto tot_pi0_mom = sqrt(pow(mom_vect[0],2) + pow(mom_vect[1],2) + pow(mom_vect[2],2) );
            //===========================================
            auto IP = pow(_geoAlgo.SqDist(shr1_bkwrd_hl,shr2_bkwrd_hl),0.5);
            auto radL_shr1 = vertex_reco.Dist(shr1.ShowerStart());
            auto radL_shr2 = vertex_reco.Dist(shr2.ShowerStart());

            auto dist_temp = sqrt( pow(v.X() - vertex.at(0),2) + pow(v.Y() - vertex.at(1),2) +pow(v.Z() - vertex.at(2),2));

            // Make the backwards projection for the showers
            auto gamm1_IP = _geoAlgo.SqDist(vertex_reco, shr1_bkwrd_hl) ;
            auto gamm2_IP = _geoAlgo.SqDist(vertex_reco, shr2_bkwrd_hl) ;

            _gamma1_vtx_IP = shr1.Energy(2) < shr2.Energy(2) ? gamm1_IP : gamm2_IP ;
            _gamma2_vtx_IP = shr1.Energy(2) < shr2.Energy(2) ? gamm2_IP : gamm1_IP ;

            //std::cout<<"DISTANCE: "<<dist_temp<<", "<<_gamma1_vtx_IP <<std::endl ;

            _pi0_mass      = sqrt(2 * shr1.Energy(2) * shr2.Energy(2) *(1.-cos(oangle))); 
            _pi0_mom       = tot_pi0_mom;
            _gamma_oangle    = oangle;
            _gamma_low_E  = shr1.Energy(2) < shr2.Energy(2) ? shr1.Energy(2) : shr2.Energy(2) ;
            _gamma_high_E = shr1.Energy(2) < shr2.Energy(2) ? shr2.Energy(2) : shr1.Energy(2) ;
            _gamma_low_RL  = shr1.Energy(2) < shr2.Energy(2) ? radL_shr1 : radL_shr2 ;
            _gamma_high_RL = shr1.Energy(2) < shr2.Energy(2) ? radL_shr2 : radL_shr1 ;

            _gamma_low_dirx = shr1.Energy(2) < shr2.Energy(2) ? dir1.X() : dir2.X() ;
            _gamma_low_diry = shr1.Energy(2) < shr2.Energy(2) ? dir1.Y() : dir2.Y() ;
            _gamma_low_dirz = shr1.Energy(2) < shr2.Energy(2) ? dir1.Z() : dir2.Z() ;
            _gamma_high_dirx = shr1.Energy(2) < shr2.Energy(2) ? dir2.X() : dir1.X() ;
            _gamma_high_diry = shr1.Energy(2) < shr2.Energy(2) ? dir2.Y() : dir1.Y() ;
            _gamma_high_dirz = shr1.Energy(2) < shr2.Energy(2) ? dir2.Z() : dir1.Z() ;

            _gamma_low_startx = shr1.Energy(2) < shr2.Energy(2) ? st1.X() : st2.X() ;
            _gamma_low_starty = shr1.Energy(2) < shr2.Energy(2) ? st1.Y() : st2.Y() ;
            _gamma_low_startz = shr1.Energy(2) < shr2.Energy(2) ? st1.Z() : st2.Z() ;
            _gamma_high_startx = shr1.Energy(2) < shr2.Energy(2) ? st2.X() : st1.X() ;
            _gamma_high_starty = shr1.Energy(2) < shr2.Energy(2) ? st2.Y() : st1.Y() ;
            _gamma_high_startz = shr1.Energy(2) < shr2.Energy(2) ? st2.Z() : st1.Z() ;

            _gamma_IP = IP;
            //_gamma1_vtx_IP = sqrt( pow( v.X() - st1.X(),2) + pow( v.Y() - st1.Y(),2) + pow( v.Z() - st1.Z(),2)); 
            //_gamma2_vtx_IP = sqrt( pow( v.X() - st2.X(),2) + pow( v.Y() - st2.Y(),2) + pow( v.Z() - st2.Z(),2)); 

            bool shower1islow = shr1.Energy(2) < shr2.Energy(2) ? true : false ;



            if ( s1 == reco_g1_id ){
              if( shower1islow){
                   auto mcs = ev_mcshr->at(mc_g1_id);
                   _mc_low_dirx = mcs.DetProfile().Px();  
                   _mc_low_diry = mcs.DetProfile().Py();  
                   _mc_low_dirz = mcs.DetProfile().Pz();  

                   //_mc_low_e = mcs.DetProfile().E() ;

                   auto mcx = mcs.DetProfile().X() ;
                   auto mcy = mcs.DetProfile().Y() ;
                   auto mcz = mcs.DetProfile().Z() ;
                   auto mct = mcs.DetProfile().T() ;
                   auto vtxtick = (mct/ 1000.) * 2.; 
                   auto vtxtimecm = vtxtick * _time2cm; 
                   auto sce_corr = _SCE->GetPosOffsets(mcx,mcy,mcz);

                   _mc_low_startx = mcx + vtxtimecm + 0.7 - sce_corr.at(0);
                   _mc_low_starty = mcy + sce_corr.at(1);
                   _mc_low_startz = mcz + sce_corr.at(2);

                  _res_low = dot_g1 ;
               }
               else{
                   auto mcs = ev_mcshr->at(mc_g1_id);
                   _mc_high_dirx = mcs.DetProfile().Px();  
                   _mc_high_diry = mcs.DetProfile().Py();  
                   _mc_high_dirz = mcs.DetProfile().Pz();  

                   //_mc_high_startx = mcs.DetProfile().X();  
                   //_mc_high_starty = mcs.DetProfile().Y();  
                   //_mc_high_startz = mcs.DetProfile().Z();  

                   auto mcx = mcs.DetProfile().X() ;
                   auto mcy = mcs.DetProfile().Y() ;
                   auto mcz = mcs.DetProfile().Z() ;
                   auto mct = mcs.DetProfile().T() ;
                   auto vtxtick = (mct/ 1000.) * 2.; 
                   auto vtxtimecm = vtxtick * _time2cm; 
                   auto sce_corr = _SCE->GetPosOffsets(mcx,mcy,mcz);

                   _mc_high_startx = mcx + vtxtimecm + 0.7 - sce_corr.at(0);
                   _mc_high_starty = mcy + sce_corr.at(1);
                   _mc_high_startz = mcz + sce_corr.at(2);

                  _res_high = dot_g1 ;
               }
             }
            
             if ( s1 == reco_g2_id ){
              if( shower1islow){
                   auto mcs = ev_mcshr->at(mc_g2_id);
                   _mc_low_dirx = mcs.DetProfile().Px();  
                   _mc_low_diry = mcs.DetProfile().Py();  
                   _mc_low_dirz = mcs.DetProfile().Pz();  

                   _mc_low_startx = mcs.DetProfile().X();  
                   _mc_low_starty = mcs.DetProfile().Y();  
                   _mc_low_startz = mcs.DetProfile().Z();  

                  _res_low = dot_g2 ;
               }
               else{
                   auto mcs = ev_mcshr->at(mc_g2_id);
                   _mc_high_dirx = mcs.DetProfile().Px();  
                   _mc_high_diry = mcs.DetProfile().Py();  
                   _mc_high_dirz = mcs.DetProfile().Pz();  

                   _mc_high_startx = mcs.DetProfile().X();  
                   _mc_high_starty = mcs.DetProfile().Y();  
                   _mc_high_startz = mcs.DetProfile().Z();  
                  _res_high = dot_g2 ;
               }
             }

           if ( s2 == reco_g1_id ){
              if( !shower1islow ){
                   auto mcs = ev_mcshr->at(mc_g1_id);
                   _mc_low_dirx = mcs.DetProfile().Px();  
                   _mc_low_diry = mcs.DetProfile().Py();  
                   _mc_low_dirz = mcs.DetProfile().Pz();  

                   _mc_low_startx = mcs.DetProfile().X();  
                   _mc_low_starty = mcs.DetProfile().Y();  
                   _mc_low_startz = mcs.DetProfile().Z();  
                  _res_low = dot_g1 ;
               }
               else{
                   auto mcs = ev_mcshr->at(mc_g1_id);
                   _mc_high_dirx = mcs.DetProfile().Px();  
                   _mc_high_diry = mcs.DetProfile().Py();  
                   _mc_high_dirz = mcs.DetProfile().Pz();  

                   _mc_high_startx = mcs.DetProfile().X();  
                   _mc_high_starty = mcs.DetProfile().Y();  
                   _mc_high_startz = mcs.DetProfile().Z();  
                  _res_high = dot_g1 ;
               }
             }
            
             if ( s2 == reco_g2_id ){
              if( !shower1islow){
                   auto mcs = ev_mcshr->at(mc_g2_id);
                   _mc_low_dirx = mcs.DetProfile().Px();  
                   _mc_low_diry = mcs.DetProfile().Py();  
                   _mc_low_dirz = mcs.DetProfile().Pz();  

                   _mc_low_startx = mcs.DetProfile().X();  
                   _mc_low_starty = mcs.DetProfile().Y();  
                   _mc_low_startz = mcs.DetProfile().Z();  
                  _res_low = dot_g2 ;
               }
               else{
                   auto mcs = ev_mcshr->at(mc_g2_id);
                   _mc_high_dirx = mcs.DetProfile().Px();  
                   _mc_high_diry = mcs.DetProfile().Py();  
                   _mc_high_dirz = mcs.DetProfile().Pz();  

                   _mc_high_startx = mcs.DetProfile().X();  
                   _mc_high_starty = mcs.DetProfile().Y();  
                   _mc_high_startz = mcs.DetProfile().Z();  
                  _res_high = dot_g2 ;
               }
             }

            auto st_high_res = sqrt ( pow(_mc_high_startx - _gamma_high_startx,2) +pow(_mc_high_starty - _gamma_high_starty,2)+ pow(_mc_high_startz - _gamma_high_startz,2) ) ;
            auto st_low_res = sqrt ( pow(_mc_low_startx - _gamma_low_startx,2) +pow(_mc_low_starty - _gamma_low_starty,2)+ pow(_mc_low_startz - _gamma_low_startz,2) ) ;
            

            if ( _gamma_IP  > 4 ){
              std::cout<<"EVENT !! "<<_event<<", IP: "<<_gamma_IP<<", "<<_gamma1_vtx_IP<<", "<<_gamma2_vtx_IP<<std::endl ;
              std::cout<<"Res: "<<_res_low<<", "<<_res_high<<", start low high : "<<st_low_res<<", "<<st_high_res<<std::endl ;

              //auto geomH = ::larutil::GeometryHelper::GetME();
              //auto high_mc_v = { _mc_high_startx, _mc_high_starty, _mc_high_startz };
              //auto low_mc_v = { _mc_low_startx, _mc_low_starty, _mc_low_startz };
              //auto high_gamma_v = { _gamma_high_startx, _gamma_high_starty, _gamma_high_startz };
              //auto low_gamma_v = { _gamma_low_startx, _gamma_low_starty, _gamma_low_startz };

              //auto high_2d_mc = geomH->Point_3Dto2D(high_mc_v,2);
              //auto low_2d_mc = geomH->Point_3Dto2D(low_mc_v,2);
              //auto high_2d_gamma = geomH->Point_3Dto2D(high_gamma_v,2);
              //auto low_2d_gamma = geomH->Point_3Dto2D(low_gamma_v,2);

              //auto w = geomH->WireToCm() ;
              //auto t = geomH->TimeToCm() ;

              //std::cout<<" High MC (w,t): "<< high_2d_mc.w/w <<", "<<high_2d_mc.t/t+800<<", Low MC : "<<low_2d_mc.w/w<<", "<<low_2d_mc.t/t +800 <<std::endl;
              //std::cout<<" High Reco (w,t): "<< high_2d_gamma.w/w <<", "<<high_2d_gamma.t/t+800<<", Low Reco: "<<low_2d_gamma.w/w<<", "<<low_2d_gamma.t/t+800 <<std::endl;


             }
            else
              std::cout<<"Normal EVENT !! "<<_event<<", IP: "<<_gamma_IP<<", "<<_gamma1_vtx_IP<<", "<<_gamma2_vtx_IP<<std::endl ;

             //////////////////////////////////////
             float mc_clus_e = 0.;
             closest_mcs_id = -1 ;
             float closest_e = 1e9 ;

             for (size_t j = 0; j < ass_showerreco_v.at(s1).size(); j++ ){

               auto clus_id = ass_showerreco_v.at(s1).at(j);
               auto iclus = ev_clus->at(clus_id);

               int plane = iclus.Plane().Plane ;
               if ( plane != 2 ) continue;

               pur_ctr_v.clear();
               pur_ctr_v.resize(ass_mcclus_v.size()+1,0) ;

               max_hits = -1;
               max_cid = -1 ;

               // Loop through all hits associared to the cluster 
               for ( int k = 0; k < ass_imageclus_v.at(clus_id).size(); k++ ){

                 auto hid = ass_imageclus_v.at(clus_id).at(k) ;
                 auto h = ev_hit->at(hid);

                 if ( _mc_hit_map.find(hid) != _mc_hit_map.end() ){

                   auto mcclus_id = _mc_hit_map[hid] ;
                   pur_ctr_v[mcclus_id]++ ;

                   if( pur_ctr_v[ mcclus_id] > max_hits ){
                     max_hits = pur_ctr_v[mcclus_id];
                     max_cid = mcclus_id ;
                   }
                 }   
                 else {
                   auto mcclus_id = ass_mcclus_v.size() ;
                   pur_ctr_v[mcclus_id]++ ; 
                   if( pur_ctr_v[mcclus_id] > max_hits ){
                     max_hits = pur_ctr_v[mcclus_id];
                     max_cid = mcclus_id ; 
                   }
                 }
               }

                if ( max_cid != ass_mcclus_v.size() && max_cid != -1 ){

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

                  // Find mcs this cluster belongs to in order to store the true shower energy
                  for ( int i = 0; i < ev_mcshr->size(); i++ ) { 

                    auto s = ev_mcshr->at(i) ;
                    if(s.Origin() != 1 || s.MotherPdgCode() != 111 ) continue;
    
                    auto e = fabs(mc_clus_e - s.DetProfile().E()) ;

                    if ( e < closest_e ){
                      closest_e  = e;
                      closest_mcs_id = i ; 
                    }   
                  }   


                  auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
                  auto tot_reco_hits = ass_imageclus_v[clus_id].size();

                  auto mcclus2 = ev_mcc->at(max_cid) ;
                  auto shr2_origin = mcclus2.Width() ;
                  auto shr2_type = mcclus2.StartOpeningAngle() ;

                  if ( shr2_origin != 2 && shr1_origin != 2 )
                    _pi0_origin = true;
                  if ( shr2_type == 1 || shr1_type == 1 )
                    _pi0_type = true;
    
                  if ( shower1islow ){ 
                    _gamma_low_purity   = _gamma_purity ;
                    _gamma_low_complete = _gamma_complete ;
                    _gamma_high_purity = float(max_hits) / tot_reco_hits ;
                    _gamma_high_complete = float(max_hits) / tot_mc_hits ;

                    //if ( closest_mcs_id != -1 ){

                    //  _mc_low_dirx = _mc_dirx ;
                    //  _mc_low_diry = _mc_diry ;
                    //  _mc_low_dirz = _mc_dirz ;
                    //  _mc_high_dirx = ev_mcshr->at(closest_mcs_id).DetProfile().Px() ;
                    //  _mc_high_diry = ev_mcshr->at(closest_mcs_id).DetProfile().Py() ;
                    //  _mc_high_dirz = ev_mcshr->at(closest_mcs_id).DetProfile().Pz() ;
                    //}

                  }
                  else{
                    _gamma_high_purity   = _gamma_purity ;
                    _gamma_high_complete = _gamma_complete ;
                    _gamma_low_purity   = float(max_hits) / tot_reco_hits ;
                    _gamma_low_complete = float(max_hits) / tot_mc_hits ;

                    //if ( closest_mcs_id != -1 ){
                    //  _mc_high_dirx = _mc_dirx ;
                    //  _mc_high_diry = _mc_diry ;
                    //  _mc_high_dirz = _mc_dirz ;
                    //  _mc_low_dirx = ev_mcshr->at(closest_mcs_id).DetProfile().Px() ;
                    //  _mc_low_diry = ev_mcshr->at(closest_mcs_id).DetProfile().Py() ;
                    //  _mc_low_dirz = ev_mcshr->at(closest_mcs_id).DetProfile().Pz() ;
                    //}

                  }

                //std::cout<<"res low and high : "<<_res_low<<", "<<_res_high<<std::endl ;
                }
             }
             ////////////////////////////////////

             mc_clus_e = 0.;
             closest_mcs_id = -1 ;
             closest_e = 1e9 ;

             for (size_t j = 0; j < ass_showerreco_v.at(s2).size(); j++ ){

               auto clus_id = ass_showerreco_v.at(s2).at(j);
               auto iclus = ev_clus->at(clus_id);

               int plane = iclus.Plane().Plane ;
               if ( plane != 2 ) continue;

               pur_ctr_v.clear();
               pur_ctr_v.resize(ass_mcclus_v.size()+1,0) ;

               max_hits = -1;
               max_cid = -1 ;

               // Loop through all hits associared to the cluster 
               for ( int k = 0; k < ass_imageclus_v.at(clus_id).size(); k++ ){

                 auto hid = ass_imageclus_v.at(clus_id).at(k) ;
                 auto h = ev_hit->at(hid);

                 if ( _mc_hit_map.find(hid) != _mc_hit_map.end() ){

                   auto mcclus_id = _mc_hit_map[hid] ;
                   pur_ctr_v[mcclus_id]++ ;

                   if( pur_ctr_v[ mcclus_id] > max_hits ){
                     max_hits = pur_ctr_v[mcclus_id];
                     max_cid = mcclus_id ;
                   }
                 }   
                 else {
                   auto mcclus_id = ass_mcclus_v.size() ;
                   pur_ctr_v[mcclus_id]++ ; 
                   if( pur_ctr_v[mcclus_id] > max_hits ){
                     max_hits = pur_ctr_v[mcclus_id];
                     max_cid = mcclus_id ; 
                   }
                 }
               }

                if ( max_cid != ass_mcclus_v.size() && max_cid != -1 ){

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

                  // Find mcs this cluster belongs to in order to store the true shower energy
                  for ( int i = 0; i < ev_mcshr->size(); i++ ) { 

                    auto s = ev_mcshr->at(i) ;
                    if(s.Origin() != 1 || s.MotherPdgCode() != 111 ) continue;
    
                    auto e = fabs(mc_clus_e - s.DetProfile().E()) ;

                    if ( e < closest_e ){
                      closest_e  = e;
                      closest_mcs_id = i ; 
                    }   
                  }   

                  auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
                  auto tot_reco_hits = ass_imageclus_v[clus_id].size();
                 if ( !shower1islow ){ 
                    _gamma_low_purity   = float(max_hits) / tot_reco_hits ;
                    _gamma_low_complete = float(max_hits) / tot_mc_hits ;
                  }
                  else{
                    _gamma_high_purity   = float(max_hits) / tot_reco_hits ;
                    _gamma_high_complete = float(max_hits) / tot_mc_hits ;
                  }
                }
             }
             ////////////////////////////////////


            _gamma_tree->Fill() ;
        }// shower ID 2 

 
      }// shower ID 1 

      //std::cout<<"JustifyPi0Cuts - Found a candidate! "<<std::endl ;
      
    return true;
  }

  bool JustifyPi0Cuts::finalize() {

    if(_fout) { 
      _fout->cd(); 
      _compare_tree->Write(); 
      _gamma_tree->Write(); 
      _one_gamma_tree->Write(); 
      _tree->Write();
    }

    std::cout<<"Signals: "<<std::endl ;
    std::cout<<"Total CCpi0 : "<<_n_cc1pi0<<std::endl; 

    // Note that cc other includes secondary pi0s.
    std::cout<<"\nBackgrounds: "<<std::endl;
    std::cout<<"1) Cosmic : "<<_n_cosmic<< std::endl;
    std::cout<<"2) CC 1pi0 : "<<_n_cc1pi0<<std::endl;
    std::cout<<"3) CC 0pi0 : "<<_n_cc0pi0<<std::endl;
    std::cout<<"4) NC 1pi0 : "<<_n_nc1pi0<<std::endl;
    std::cout<<"5) NC 0pi0 : "<<_n_nc0pi0<<std::endl;
    std::cout<<"6) Other   : "<<_n_other<<std::endl; 

    std::cout<<"Total accounted backgrounds: "<< _n_other + _n_cosmic + _n_nc1pi0 + _n_nc0pi0 + _n_cc0pi0 <<std::endl ;

    return true;
  }

}
#endif
