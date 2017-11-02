#ifndef LARLITE_JUSTIFYPI0CUTS_CXX
#define LARLITE_JUSTIFYPI0CUTS_CXX

#include "JustifyPi0Cuts.h"
#include "DataFormat/shower.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/cluster.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"

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
    }

    if( !_one_gamma_tree ){
      _one_gamma_tree = new TTree("one_gamma_tree","");
      _one_gamma_tree->Branch("_event",&_event,"event/I");
      _one_gamma_tree->Branch("_gamma_E",&_gamma_E,"gamma_E/F");
      _one_gamma_tree->Branch("_gamma_RL",&_gamma_RL,"gamma_RL/F");
      _one_gamma_tree->Branch("_gamma_vtx_IP",&_gamma_vtx_IP,"gamma_vtx_IP/F");
      _one_gamma_tree->Branch("_gamma_matched",&_gamma_matched,"gamma_matched/B");
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

    //std::cout<<"\nEvent : "<<_event <<std::endl;

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
      _mu_mom  = t.VertexMomentum() ;
      _mu_len  = t.Length(0); // Calculates the length from point 0 to end
      _mu_angle = t.Theta() ;
      _mu_phi = t.Phi() ; 

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
      auto vtxtimecm = vtxtick * _time2cm; // time in cm :

      //Map of lengths -> track id
      std::multimap<float,int> trk_map ;

      // Grab the origin of the track and assess backgrounds properly
      std::multimap<float,int> mctrk_map ;
      auto tag_st = t.Vertex() ;
      auto tag_end = t.End() ;

      for ( size_t ti = 0; ti < ev_mctrk->size(); ti++ ) { 

        auto mc_vtx = ev_mctrk->at(ti).Start() ;
        auto mc_end = ev_mctrk->at(ti).End() ;
      
        float dist_st = sqrt(  pow(mc_vtx.X() - _mu_startx,2) + 
                               pow(mc_vtx.Y() - _mu_starty,2) + 
                               pow(mc_vtx.Z() - _mu_startz,2) );  

        float dist_end = sqrt( pow(mc_vtx.X() - _mu_endx,2) + 
                               pow(mc_vtx.Y() - _mu_endx,2) + 
                               pow(mc_vtx.Z() - _mu_endx,2) );  

         
         if ( dist_st < 25 || dist_end < 25){
            float length = sqrt( pow(mc_end.X() - mc_vtx.X(),2) + 
                             pow(mc_end.Y() - mc_vtx.Y(),2) + 
                             pow(mc_end.Z() - mc_vtx.Z(),2) );  

            mctrk_map.emplace(1./length,ti);
         }   
       }   
       
       int mc_max_it = -1;
       float mc_max_dot = -1.;

       if( mctrk_map.size() ) { 

        auto tag_st = t.VertexDirection();     
        auto tag_norm = sqrt( pow(tag_st.Px(),2) + pow(tag_st.Py(),2) + pow(tag_st.Pz(),2)); 

        for( auto & ti : mctrk_map ){
              
          auto mc = ev_mctrk->at(ti.second);
          auto mc_st = mc.Start();
          auto mc_norm = sqrt( pow(mc_st.Px(),2) + pow(mc_st.Py(),2) + pow(mc_st.Pz(),2) );
          
          auto dot = (tag_st.Px() * mc_st.Px() + tag_st.Py() * mc_st.Py() + tag_st.Pz() * mc_st.Pz())/tag_norm / mc_norm ;

          if ( fabs(dot) > mc_max_dot ){
               mc_max_dot = dot;
               mc_max_it = ti.second ;
          }
        }
      }   
      // If no true tracks aligned with reco track, mark it as cosmic 
      else {
        _n_cosmic++;
        _bkgd_id = 1 ;
      }

      if( _bkgd_id == -1 ){

        bool infv = true;
        if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
          infv = false;

        auto parts = ev_mctruth->at(0).GetParticles();
        int n_pi0 = 0;

        if( ev_mctrk->at(mc_max_it).Origin() == 2 ){
          _n_cosmic++;
          _bkgd_id = 1; 
        }
        for ( auto const & p : parts ){
          if( p.StatusCode() == 1 && p.PdgCode() == 111 )
            n_pi0 ++;
        }   

        if( nu.CCNC() == 0 && n_pi0 == 1 && infv ) {
          _bkgd_id = 2;
          _n_cc1pi0 ++; 
        }
        else if( nu.CCNC() == 0 && n_pi0 == 0 ) {
          _bkgd_id = 3;
          _n_cc0pi0++;
        }
        else if( nu.CCNC() == 1 && n_pi0 == 1 ) {
          _bkgd_id = 4;
          _n_nc1pi0 ++; 
        }
        else if( nu.CCNC() == 1 && n_pi0 == 0 ) {
          _bkgd_id = 5;
          _n_nc0pi0++;
        }
        else {
          _bkgd_id = 6;
          _n_other ++;   

        }
      }

    auto ev_s = storage->get_data<event_shower>("showerreco");
    _nshrs = ev_s->size() ;
   
    _tree->Fill();

    if( !ev_s || !ev_s->size() || ev_s->size() < 1 ){
      std::cout<<"Not enough reco'd showers..." <<ev_s->size()<<std::endl;
      return false;
    } 

    auto ev_mcshr = storage->get_data<event_mcshower>("mcreco");

    // Now deal with shower comparisons
    std::vector<int> shr_ids;

    bool found_pi0 = false;

    for ( int si = 0; si < ev_mcshr->size(); si++){

      auto s = ev_mcshr->at(si);
      auto st = s.Start();
      auto dist = sqrt( pow(st.X() - xyz[0],2) + pow(st.Y() - xyz[1],2) + pow(st.Z() - xyz[2],2) );

      if ( s.DetProfile().E() > 0 && s.MotherPdgCode() == 111 && s.Origin() == 1 ) found_pi0 = true;

      if ( dist < 0.00001 && found_pi0 && nu.Nu().PdgCode() == 14 && nu.CCNC() == 0 ){
        shr_ids.emplace_back(si) ;
        _event_type = 0;
      }
    }

    if ( found_pi0 && !shr_ids.size() ){
      for ( int si = 0; si < ev_mcshr->size(); si++){

        auto s = ev_mcshr->at(si);
        auto st = s.Start();
        auto dist = sqrt( pow(st.X() - xyz[0],2) + pow(st.Y() - xyz[1],2) + pow(st.Z() - xyz[2],2) );

        if ( s.DetProfile().E() > 0 && s.MotherPdgCode() == 111 ){
          shr_ids.emplace_back(si) ;
	  _event_type = 1;
	  std::cout<<"FOUND A 1 event!" <<std::endl ;
	  }
	}
    }
      
    if ( _event_type == -1 ) _event_type = 2;

    std::multimap<float,std::pair<int,int>> mc_reco_map ;    

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

        if ( fabs(dot) > 1 ) std::cout<<"DOT ! " <<dot <<std::endl ;

        mc_reco_map.emplace(1./dot,std::make_pair(mc_id,reco_id)) ;
      }
    } 

    int reco_g1_id = -1, reco_g2_id = -1;
    int mc_g1_id = -1, mc_g2_id = -1;
    float dot_g1 = -10, dot_g2 = -10;

    for ( auto const & m : mc_reco_map ){

      auto score = 1./m.first ;
      std::cout<<"Some scores are..." <<score<<std::endl;
    
      if ( score < 0.9 ) break; 

      if ( reco_g1_id == -1 ){
        mc_g1_id   = m.second.first ;
        reco_g1_id = m.second.second ;
        dot_g1 = score ;
      }
      else{
        if ( m.second.first == mc_g1_id || m.second.second == reco_g1_id ) continue;
        
        mc_g2_id = m.second.first;       
        reco_g2_id = m.second.second ;       
        dot_g2 = score ;
        break;
      }
    }

    for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

       clearCompare();

       auto const& shr1 = ev_s->at(s1);

       _gamma_E  = shr1.Energy()  ;
       ::geoalgo::Point_t vertex_reco(v.X(),v.Y(),v.Z());
        _gamma_RL = vertex_reco.Dist(shr1.ShowerStart());

       geoalgo::Vector_t rev_shr(-1.*shr1.Direction()) ;
       auto shr_bkwrd_hl = ::geoalgo::HalfLine_t(shr1.ShowerStart(),rev_shr);
       _gamma_vtx_IP = _geoAlgo.SqDist(vertex_reco, shr_bkwrd_hl) ;

       _reco_E = _gamma_E ;
       _reco_startx = shr1.ShowerStart().X();  
       _reco_starty = shr1.ShowerStart().Y();  
       _reco_startz = shr1.ShowerStart().Z();  
       _reco_dirx = shr1.Direction().Px();  
       _reco_diry = shr1.Direction().Py();  
       _reco_dirz = shr1.Direction().Pz();  

       if ( s1 == reco_g1_id ){
          _gamma_high_matched = true ;
          _gamma_matched = true; 

	 auto mcs = ev_mcshr->at(mc_g1_id);

         auto st_x = mcs.DetProfile().X();
         auto st_y = mcs.DetProfile().Y();
         auto st_z = mcs.DetProfile().Z();

         auto sce_corr = _SCE->GetPosOffsets(st_x,st_y,st_z);
         _mc_startx = st_x + vtxtimecm - sce_corr.at(0); 
         _mc_starty = st_y + sce_corr.at(1); 
         _mc_startz = st_z + sce_corr.at(2);

         _mc_E = mcs.DetProfile().E();  
         _mc_dirx = mcs.DetProfile().Px();  
         _mc_diry = mcs.DetProfile().Py();  
         _mc_dirz = mcs.DetProfile().Pz();  

	 _reco_dot = dot_g1 ;
	 _reco_start3D = sqrt( pow(_mc_startx - _reco_startx,2) + pow(_mc_starty - _reco_starty,2) + pow(_mc_startz - _reco_startz,2)); 

         _compare_tree->Fill();
       }
       else if ( s1 == reco_g2_id ){

	 auto mcs = ev_mcshr->at(mc_g2_id);

         auto st_x = mcs.DetProfile().X();
         auto st_y = mcs.DetProfile().Y();
         auto st_z = mcs.DetProfile().Z();

         auto sce_corr = _SCE->GetPosOffsets(st_x,st_y,st_z);
         _mc_startx = st_x + vtxtimecm - sce_corr.at(0); 
         _mc_starty = st_y + sce_corr.at(1); 
         _mc_startz = st_z + sce_corr.at(2);

         _mc_E = mcs.DetProfile().E();  
         _mc_dirx = mcs.DetProfile().Px();  
         _mc_diry = mcs.DetProfile().Py();  
         _mc_dirz = mcs.DetProfile().Pz();  
	 
	 _reco_dot = dot_g2 ;
	 _reco_start3D = sqrt( pow(_mc_startx - _reco_startx,2) + pow(_mc_starty - _reco_starty,2) + pow(_mc_startz - _reco_startx,2)); 

         _compare_tree->Fill();
       }


       _one_gamma_tree->Fill();

       for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

            if (s2 <= s1) continue;

            if ( s2 == reco_g2_id ) _gamma_low_matched = true ;

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
            _gamma1_vtx_IP = _geoAlgo.SqDist(vertex_reco, shr1_bkwrd_hl) ;
            _gamma2_vtx_IP = _geoAlgo.SqDist(vertex_reco, shr2_bkwrd_hl) ;

            //std::cout<<"DISTANCE: "<<dist_temp<<", "<<_gamma1_vtx_IP <<std::endl ;

            _pi0_mass      = sqrt(2 * shr1.Energy() * shr2.Energy() *(1.-cos(oangle))); 
            _pi0_mom       = tot_pi0_mom;
            _gamma_oangle    = oangle;
            _gamma_low_E  = shr1.Energy() < shr2.Energy() ? shr1.Energy() : shr2.Energy() ;
            _gamma_high_E = shr1.Energy() < shr2.Energy() ? shr2.Energy() : shr1.Energy() ;
            _gamma_low_RL  = shr1.Energy() < shr2.Energy() ? radL_shr1 : radL_shr2 ;
            _gamma_high_RL = shr1.Energy() < shr2.Energy() ? radL_shr2 : radL_shr1 ;
            _gamma_IP = IP;
            _gamma1_vtx_IP = sqrt( pow( v.X() - st1.X(),2) + pow( v.Y() - st1.Y(),2) + pow( v.Z() - st1.Z(),2)); 
            _gamma2_vtx_IP = sqrt( pow( v.X() - st2.X(),2) + pow( v.Y() - st2.Y(),2) + pow( v.Z() - st2.Z(),2)); 

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
