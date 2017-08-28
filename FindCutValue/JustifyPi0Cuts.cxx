#ifndef LARLITE_JUSTIFYPI0CUTS_CXX
#define LARLITE_JUSTIFYPI0CUTS_CXX

#include "JustifyPi0Cuts.h"
#include "DataFormat/shower.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/cluster.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"

namespace larlite {

  bool JustifyPi0Cuts::initialize() {

    _n_other = 0;    // 0 
    _n_cosmic = 0;   // 1
    _n_cc1pi0 = 0;   // 2 
    _n_cc0pi0 = 0;   // 3
    _n_nc1pi0 = 0;   // 4 
    _n_nc0pi0 = 0;   // 5

    
    if( !_pi0_selection ){
      _pi0_selection = new TTree("pi0_selection","");
      _pi0_selection->Branch("_event",&_event,"event/I");
      _pi0_selection->Branch("_pi0_mass",&_pi0_mass,"pi0_mass/F");
      _pi0_selection->Branch("_pi0_mom",&_pi0_mom,"pi0_mom/F");
      _pi0_selection->Branch("_pi0_oangle",&_pi0_oangle,"pi0_oangle/F");
      _pi0_selection->Branch("_pi0_low_shrE",&_pi0_low_shrE,"pi0_low_shrE/F");
      _pi0_selection->Branch("_pi0_high_shrE",&_pi0_high_shrE,"pi0_high_shrE/F");
      _pi0_selection->Branch("_pi0_low_radL",&_pi0_low_radL,"pi0_low_radL/F");
      _pi0_selection->Branch("_pi0_high_radL",&_pi0_high_radL,"pi0_high_radL/F");
      _pi0_selection->Branch("_pi0_IP",&_pi0_IP,"pi0_IP/F");
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
    }

    _event = -1;

    return true;
  }

  void JustifyPi0Cuts::clear(){
  
    _pi0_mass      = -10;
    _pi0_mom       = -10;
    _pi0_oangle    = -10;
    _pi0_low_shrE  = -10;
    _pi0_high_shrE = -10;
    _pi0_low_radL  = -10;
    _pi0_high_radL = -10;
    _pi0_IP = -10;
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
  }
  
  bool JustifyPi0Cuts::analyze(storage_manager* storage) {

    _event++;

    auto ev_s = storage->get_data<event_shower>("showerreco");

    if( !ev_s || !ev_s->size() || ev_s->size() < 2 ){

      std::cout<<"Not enough reco'd showers..." <<ev_s->size()<<std::endl;
      return false;
     }

    std::cout<<"\nEvent : "<<_event <<std::endl;

    clear();

    std::vector<std::pair<int,int>> candidate_pairs;
    std::vector<int> cand_ids;

    for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

        auto const& shr1 = ev_s->at(s1);

        for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

            if (s2 <= s1) continue;

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
            auto radL_shr1 = vertex.Dist(shr1.ShowerStart());
            auto radL_shr2 = vertex.Dist(shr2.ShowerStart());

            //if( oangle < 0.35 ) continue;
            //if( pow( _geoAlgo.SqDist(shr1_bkwrd_hl, shr2_bkwrd_hl), 0.5 ) > 4.) continue; 
            //if( radL_shr1 > 62 || radL_shr2 > 62 ) continue;

            _pi0_mass      = sqrt(2 * shr1.Energy() * shr2.Energy() *(1.-cos(oangle))); 
            _pi0_mom       = tot_pi0_mom;
            _pi0_oangle    = oangle;
            _pi0_low_shrE  = shr1.Energy() < shr2.Energy() ? shr1.Energy() : shr2.Energy() ;
            _pi0_high_shrE = shr1.Energy() < shr2.Energy() ? shr2.Energy() : shr1.Energy() ;
            _pi0_low_radL  = shr1.Energy() < shr2.Energy() ? radL_shr1 : radL_shr2 ;
            _pi0_high_radL = shr1.Energy() < shr2.Energy() ? radL_shr2 : radL_shr1 ;
            _pi0_IP = IP;
	    _pi0_selection->Fill() ;
        }// shower ID 2 
      }// shower ID 1 

      std::cout<<"JustifyPi0Cuts - Found a candidate! "<<std::endl ;

      auto ev_v = storage->get_data<event_vertex>("numuCC_vertex");
      auto ev_t = storage->get_data<event_track>("numuCC_track");

      if( !ev_t || !ev_t->size() ){
        std::cout<<"No tagged track; what??" <<std::endl;
        return false;
      }

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

      //Map of lengths -> track id
      std::multimap<float,int> trk_map ;

      // Grab the origin of the track and assess backgrounds properly
      std::multimap<float,int> mctrk_map ;
      auto tag_st = t.Vertex() ;
      auto tag_end = t.End() ;
      float mc_min_dist = 1e9;

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
            mc_min_dist = dist_st < dist_end ? dist_st : dist_end ; 
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
         std::cout<<"\nEvent is : "<<_event <<", mult: "<<trk_map.size()<<", "<<storage->event_id()<<", "<<storage->subrun_id()<<std::endl ;
        _n_cosmic++;
        _bkgd_id = 1 ;
        _tree->Fill();

       return false;
      }

      auto mc_vtx = ev_mctrk->at(mc_max_it).Start() ;
      auto mc_end = ev_mctrk->at(mc_max_it).End() ;
      bool infv = true;

      if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
        infv = false;

      ///////////////////////////////////////////////////////////////////////////////////////////////
      /// Now count the number of backgrounds and signals
      ///////////////////////////////////////////////////////////////////////////////////////////////
      auto parts = ev_mctruth->at(0).GetParticles();
      int n_pi0 = 0;

      if( ev_mctrk->at(mc_max_it).Origin() == 2 ){
        _n_cosmic++;
        _bkgd_id = 1; 
      }

      if( _bkgd_id == -1 ){
        for ( auto const & p : parts ){
          if( p.StatusCode() == 1 && p.PdgCode() == 111 )
            n_pi0 ++;
        }   

        if( nu.CCNC() == 0 && n_pi0 == 1 && infv ) {
          _bkgd_id = 2;
          _n_cc1pi0 ++; 
          _event_list.emplace_back(_event);
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

      _tree->Fill();
      
      _event_list.emplace_back(_event); 

    return true;
  }

  bool JustifyPi0Cuts::finalize() {

    if(_fout) { 
      _fout->cd(); 
      _pi0_selection->Write(); 
      _tree->Write();
    }

    std::cout<<"Signals: "<<std::endl ;
    std::cout<<"Total CCpi0 : "<<_n_cc1pi0<<"/"<<_event_list.size()<<std::endl; 

    // Note that cc other includes secondary pi0s.
    std::cout<<"\nBackgrounds: "<<std::endl;
    std::cout<<"1) Cosmic : "<<_n_cosmic<< std::endl;
    std::cout<<"2) CC 1pi0 : "<<_n_cc1pi0<<std::endl;
    std::cout<<"3) CC 0pi0 : "<<_n_cc0pi0<<std::endl;
    std::cout<<"4) NC 1pi0 : "<<_n_nc1pi0<<std::endl;
    std::cout<<"5) NC 0pi0 : "<<_n_nc0pi0<<std::endl;
    std::cout<<"6) Other   : "<<_n_other<<std::endl; 

    std::cout<<"Total accounted backgrounds: "<< _n_other + _n_cosmic + _n_nc1pi0 + _n_nc0pi0 + _n_cc0pi0 <<std::endl ;


   //std::cout<<"\n\n****** "<<_event_list.size()<<" events found by 2 Shower Pi0Reco Module! ******"<<std::endl; 
   //for ( auto const & e : _event_list )
   //  std::cout<<e <<", "; 

   //std::cout<<"\n\n\n";
  
    return true;
  }

}
#endif
