#ifndef LARLITE_BACKGROUNDALL_CXX
#define LARLITE_BACKGROUNDALL_CXX

#include "BackgroundAll.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/track.h"
#include "DataFormat/shower.h"
#include "DataFormat/vertex.h"

#include "LArUtil/GeometryHelper.h"

namespace larlite {

  bool BackgroundAll::initialize() {    

    _event = -1; 

    _n_other = 0;    // 0 
    _n_cosmic = 0;   // 1
    _n_cc1pi0 = 0;   // 2 
    _n_cc0pi0 = 0;   // 3
    _n_nc1pi0 = 0;   // 4 
    _n_nc0pi0 = 0;   // 5

    if ( !_tree ){
      _tree = new TTree("tree","");
      _tree->Branch("event",&_event,"event/I");
      _tree->Branch("bkgd_id",&_bkgd_id,"bkgd_id/I");

      _tree->Branch("vtx_x",&_vtx_x,"vtx_x/F");
      _tree->Branch("vtx_y",&_vtx_y,"vtx_y/F");
      _tree->Branch("vtx_z",&_vtx_z,"vtx_z/F");

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

      _tree->Branch("pi0_mass",&_pi0_mass,"pi0_mass/F");
      _tree->Branch("pi0_oangle",&_pi0_oangle,"pi0_oangle/F");
      _tree->Branch("pi0_mom",&_pi0_mom,"pi0_mom/F");
      _tree->Branch("pi0_low_shrE",&_pi0_low_shrE,"pi0_low_shrE/F");
      _tree->Branch("pi0_high_shrE",&_pi0_high_shrE,"pi0_high_shrE/F");
      _tree->Branch("pi0_low_radL",&_pi0_low_radL,"pi0_low_radL/F");
      _tree->Branch("pi0_high_radL",&_pi0_high_radL,"pi0_high_radL/F");

      _tree->Branch("gamma_E",&_gamma_E,"gamma_E/F");
      _tree->Branch("gamma_RL",&_gamma_RL,"gamma_RL/F");

   }

    return true;
  }

  void BackgroundAll::clear(){

    _bkgd_id = -1 ;
    _mult    = 0;
    _vtx_x   = -999;
    _vtx_y   = -999;
    _vtx_z   = -999;

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

    _pi0_mass = -999;
    _pi0_oangle = -999;
    _pi0_mom = -999;
    _pi0_low_shrE = -999;
    _pi0_high_shrE = -999;
    _pi0_low_radL = -999;
    _pi0_high_radL = -999;

    _gamma_E = -999;
    _gamma_RL = -999;

  }
  
  bool BackgroundAll::analyze(storage_manager* storage) {

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

    auto trk = ev_tagged_trk->at(0) ;

    // Fill track information
    _mu_phi = trk.Phi();
    _mu_angle = cos(trk.Theta());
    _mu_len =   trk.Length(0); // Allulates the length from point 0 to end
    _mu_startx = trk.Vertex().X(); 
    _mu_starty = trk.Vertex().Y(); 
    _mu_startz = trk.Vertex().Z(); 
    _mu_endx = trk.End().X(); 
    _mu_endy = trk.End().Y(); 
    _mu_endz = trk.End().Z(); 

    _mu_mom  = trk.VertexMomentum() ;

    // Adjust for pandora bug
    std::vector<double> dir = { (_mu_endx - _mu_startx) / _mu_len,
                                (_mu_endy - _mu_starty) / _mu_len,
                                (_mu_endz - _mu_startz) / _mu_len };

    auto dir_start = trk.VertexDirection();
    std::vector<double> other_dir = { dir_start.X(), dir_start.Y(), dir_start.Z() };  

    float dotProd = dir.at(0) * other_dir.at(0) + dir.at(1) * other_dir.at(1) +  dir.at(2) * other_dir.at(2) ;

    if( dotProd < 0 ) { 
       TVector3 new_dir(-dir_start.X(),-dir_start.Y(),-dir_start.Z());
       _mu_angle = cos(new_dir.Theta());
       _mu_phi = new_dir.Phi();
    }   

    // Fill multiplicity info 
    auto ev_trk = storage->get_data<event_track>("pandoraNu");

    for ( auto const & t : *ev_trk ){
        auto st = t.Vertex() ;
        auto end = t.End() ;

        auto dist_st = sqrt( pow(st.X() - vtx.X(),2) + pow(st.Y() - vtx.Y(),2) + pow(st.Z() - vtx.Z(),2) );
        auto dist_end = sqrt( pow(end.X() - vtx.X(),2) + pow(end.Y() - vtx.Y(),2) + pow(end.Z() - vtx.Z(),2) );

        if (dist_st < 3 || dist_end < 3)
          _mult ++ ;
   }

   if ( _mult == 0 ){ 
     std::cout<<"Weird...Origin? "<<_event<<", "<<storage->run_id()<<", "<<storage->subrun_id()<<", "<<storage->event_id()<<std::endl ;

     std::cout<<"Track size: " <<ev_trk->size() <<std::endl;

   }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Want to be able to access the origin of the tagged muon. Thus, need to find it, and 
    // Ask for its origin.  Need to match to MCtrack to do this
    if ( _mc_sample ){

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

      std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };

      auto vtx_diff = sqrt(pow(xyz[0] - _vtx_x,2) + pow(xyz[1] - _vtx_y,2) + pow(xyz[2] - _vtx_z,2));

      //Map of lengths -> track id
      std::multimap<float,int> trk_map ;

      // Grab the origin of the track and assess backgrounds properly
      std::multimap<float,int> mctrk_map ;
      auto tag_trk = ev_tagged_trk->at(0);
      auto tag_st = tag_trk.Vertex() ;
      auto tag_end = tag_trk.End() ;
      float mc_min_dist = 1e9;

      auto geomH = ::larutil::GeometryHelper::GetME(); 

      auto reco2d = geomH->Point_3Dto2D(tag_st,2);
      auto reco2d_end = geomH->Point_3Dto2D(tag_end,2);

      for ( size_t ti = 0; ti < ev_mctrk->size(); ti++ ) { 

        auto mc_vtx = ev_mctrk->at(ti).Start() ;
        auto mc_end = ev_mctrk->at(ti).End() ;
      
        float dist_st = sqrt(  pow(mc_vtx.X() - tag_st.X(),2) + 
                               pow(mc_vtx.Y() - tag_st.Y(),2) + 
                               pow(mc_vtx.Z() - tag_st.Z(),2) );  

        float dist_end = sqrt( pow(mc_vtx.X() - tag_end.X(),2) + 
                               pow(mc_vtx.Y() - tag_end.Y(),2) + 
                               pow(mc_vtx.Z() - tag_end.Z(),2) );  

         
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

        auto tag_st = tag_trk.VertexDirection();     
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
         //std::cout<<"\nEvent is : "<<_event <<", mult: "<<trk_map.size()<<", "<<storage->event_id()<<", "<<storage->subrun_id()<<std::endl ;
        _n_cosmic++;
        _bkgd_id = 1 ;

      }

      ///////////////////////////////////////////////////////////////////////////////////////////////
      /// Now count the number of backgrounds and signals
      ///////////////////////////////////////////////////////////////////////////////////////////////
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
    }

    if ( _get_pi0_info ){

      auto ev_s = storage->get_data<event_shower>("pi0_candidate_showers");

      if( !ev_s || !ev_s->size() || ev_s->size() < 2 ){
        std::cout<<"Not enough reco'd showers..." <<std::endl;
        return false;
       }   

      auto const& shr1 = ev_s->at(0) ;
      auto const& shr2 = ev_s->at(1);

      geoalgo::Vector_t rev_shr1(-1.*shr1.Direction()) ;
      geoalgo::Vector_t rev_shr2(-1.*shr2.Direction()) ;

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
      _pi0_low_shrE  = shr1.Energy() < shr2.Energy() ? shr1.Energy() : shr2.Energy() ;
      _pi0_high_shrE = shr1.Energy() < shr2.Energy() ? shr2.Energy() : shr1.Energy() ;
      _pi0_low_radL  = shr1.Energy() < shr2.Energy() ? radL_shr1 : radL_shr2 ;
      _pi0_high_radL = shr1.Energy() < shr2.Energy() ? radL_shr2 : radL_shr1 ;

    }
 
    if ( _get_single_shower_info ){
      auto ev_s = storage->get_data<event_shower>("pi0_1gamma_candidate_showers");
      if ( !ev_s ) return false;

      _gamma_E = ev_s->at(0).Energy(2); 
      _gamma_RL = ev_s->at(0).Energy(2); 

    }
    
    _tree->Fill();    

    return true;
  }

  bool BackgroundAll::finalize() {

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

    if ( _fout ){
      _fout->cd();
      _tree->Write();
    }
  
    return true;
  }

}
#endif
