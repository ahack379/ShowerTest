#ifndef LARLITE_SINGLEGAMMABACKGROUNDCALC_CXX
#define LARLITE_SINGLEGAMMABACKGROUNDCALC_CXX

#include "SingleGammaBackgroundCalc.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/shower.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"

namespace larlite {

  bool SingleGammaBackgroundCalc::initialize() {    

    _event = -1; 

    _n_noise = 0;     // 1
    _n_cosmic = 0;    // 2
    _n_nue = 0;       // 3
    _n_antinumu = 0;  // 4
    _n_nc = 0;        // 5
    _multpi0 = 0;     // 6
    _ccpi0_outfv = 0; // 7 
    _tot_ccpi0 = 0;   // 8
    _n_gammas = 0;    // 9
    _n_ccother = 0;   // 10

    if ( !_tree ){
      _tree = new TTree("tree","");
      _tree->Branch("event",&_event,"event/I");
      _tree->Branch("bkgd_id",&_bkgd_id,"bkgd_id/I");

      _tree->Branch("vtx_x",&_vtx_x,"vtx_x/F");
      _tree->Branch("vtx_y",&_vtx_y,"vtx_y/F");
      _tree->Branch("vtx_z",&_vtx_z,"vtx_z/F");
      _tree->Branch("mu_angle",&_mu_angle,"mu_angle/F");
      _tree->Branch("mu_len",&_mu_len,"mu_len/F");

      _tree->Branch("pi0_mass",&_pi0_mass,"pi0_mass/F");
      _tree->Branch("pi0_oangle",&_pi0_oangle,"pi0_oangle/F");
      _tree->Branch("pi0_mom",&_pi0_mom,"pi0_mom/F");
      _tree->Branch("pi0_low_shrE",&_pi0_low_shrE,"pi0_low_shrE/F");
      _tree->Branch("pi0_high_shrE",&_pi0_high_shrE,"pi0_high_shrE/F");
      _tree->Branch("pi0_low_radL",&_pi0_low_radL,"pi0_low_radL/F");
      _tree->Branch("pi0_high_radL",&_pi0_high_radL,"pi0_high_radL/F");
   }

  //std::cout<<"PI0 LIST!  "<<_pi0_list.size()<<std::endl;

    return true;
  }
  
  bool SingleGammaBackgroundCalc::analyze(storage_manager* storage) {

    auto it = _map.find(storage->subrun_id());
    bool foundit = false;

    if( it != _map.end() ){
     while ( it->first == storage->subrun_id() ){  
       auto temp_event = it->second ; 
       if( temp_event == storage->event_id() )
         foundit = true;

       it++; 
       }   
      if ( !foundit)
       _map.emplace(storage->subrun_id(), storage->event_id() );

      else return false ;
     }   
    else 
       _map.emplace(storage->subrun_id(), storage->event_id() );

     _bkgd_id = -1 ;
     _event++ ;



    //std::cout<<"\nEvent is : "<<_event <<", "<<storage->event_id()<<", "<<storage->subrun_id()<<std::endl ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
      }

    auto ev_mctrk = storage->get_data<event_mctrack>("mcreco");
    if ( !ev_mctrk || !ev_mctrk->size() ) {std::cout<<"No MCTrack!" <<std::endl ; return false; }

    auto ev_vtx= storage->get_data<event_vertex>("numuCC_vertex"); 
    if(!ev_vtx || !ev_vtx->size() ) { std::cout<<"Event has no recovertex info "<<std::endl; return false; }

    auto ev_trk = storage->get_data<event_track>("pandoraNu");
    if ( !ev_trk || !ev_trk->size() ) {std::cout<<"No Track!" <<std::endl ; return false; }

    auto ev_tagged_trk = storage->get_data<event_track>("numuCC_track");
    if ( !ev_tagged_trk || !ev_tagged_trk->size() ){ std::cout<<"No Tagged Track!" <<std::endl ; return false; }

    //auto ev_s = storage->get_data<event_shower>("showerreco");
    //auto ev_s = storage->get_data<event_shower>("pi0_1gamma_candidate_showers");

    //if( !ev_s || !ev_s->size() ) {
    //  std::cout<<"Not enough reco'd showers..." <<std::endl; return false;  }   

    // std::multimap<float,float> cand_map ;

    // for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

    //   auto vtx = ev_vtx->at(0);
    //   ::geoalgo::Point_t vertex(vtx.X(),vtx.Y(),vtx.Z());

    //   auto const& shr1 = ev_s->at(s1);
    //   geoalgo::Vector_t rev_shr1(-1.*shr1.Direction()) ;
    //   auto shr1_bkwrd_hl = ::geoalgo::HalfLine_t(shr1.ShowerStart(),rev_shr1);

    //   auto dist = _geoAlgo.SqDist(vertex, shr1_bkwrd_hl) ;

    //   auto shr_st = shr1.ShowerStart(); 
    //   auto dist_to_vtx = sqrt( pow( vtx.X() - shr_st.X(),2) +
    //                            pow( vtx.Y() - shr_st.Y(),2) + 
    //                            pow( vtx.Z() - shr_st.Z(),2)); 

    //   auto radL_shr1 = vertex.Dist(shr1.ShowerStart());
    //   _pi0_high_radL = radL_shr1 ;
    //   _pi0_high_shrE = shr1.Energy() ;

    //   if( dist > 4 ){  
    //     std::cout<<"Impact Parameter: "<<dist<<std::endl ;
    //     continue;
    //   }   

    //    if( radL_shr1 > 60 ){
    //      std::cout<<"Radiation length : "<<radL_shr1<<std::endl ;
    //      continue;
    //    }

    //    if( _pi0_high_shrE < 90 ){
    //      std::cout<<"Energy: "<<_pi0_high_shrE<<std::endl ;
    //      continue;
    //    }
    //   cand_map.emplace(dist,s1);
    // }// shower ID 1 

    // if( cand_map.size() < 1 ) return false;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Want to be able to access the origin of the tagged muon. Thus, need to find it, and 
    // Ask for its origin.  Need to match to MCtrack to do this

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    auto vtx = ev_vtx->at(0); 
    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };
    _vtx_x = vtx.X(); _vtx_y = vtx.Y(); _vtx_z = vtx.Z();

    std::multimap<float,int> trk_map ;

    // Find closest + longest pandoraNu track to vertex
    for ( size_t ti = 0; ti < ev_trk->size(); ti++ ) { 

      auto t_vtx = ev_trk->at(ti).Vertex() ;
      auto t_end = ev_trk->at(ti).End() ;
      
      float dist_st = sqrt( pow(t_vtx.X() - vtxXYZ[0],2) + 
                            pow(t_vtx.Y() - vtxXYZ[1],2) + 
                            pow(t_vtx.Z() - vtxXYZ[2],2) );  

      float dist_end = sqrt( pow(t_end.X() - vtxXYZ[0],2) + 
                             pow(t_end.Y() - vtxXYZ[1],2) + 
                             pow(t_end.Z() - vtxXYZ[2],2) );  

       if ( dist_st < 3 || dist_end < 3 ){

          float len = ev_trk->at(ti).Length();
          trk_map.emplace(1./len,ti);
       }   
     }   

    int trk_it = -1; 

    if( trk_map.size() ) 
      trk_it = trk_map.begin()->second ;
    else{
      std::cout<<"Should NOT be returning here...fix it"<<std::endl ;
      return false ;
    }

    // MC, to grab the origin of the track and assess backgrounds properly
    std::multimap<float,int> mctrk_map ;
    //auto tag_trk = ev_tagged_trk->at(0);
    auto tag_trk = ev_trk->at(trk_it);

    _mu_angle      = cos(tag_trk.Theta());
    _mu_len        = tag_trk.Length(0); // Calculates the length from point 0 to end

    auto tag_st = tag_trk.Vertex() ;
    auto tag_end = tag_trk.End() ;
    float mc_min_dist = 1e9;

    for ( size_t ti = 0; ti < ev_mctrk->size(); ti++ ) { 

      auto mc_vtx = ev_mctrk->at(ti).Start() ;
      auto mc_end = ev_mctrk->at(ti).End() ;
    
      float dist_st = sqrt(  pow(mc_vtx.X() - tag_st.X(),2) + 
                             pow(mc_vtx.Y() - tag_st.Y(),2) + 
                             pow(mc_vtx.Z() - tag_st.Z(),2) );  

      float dist_end = sqrt( pow(mc_vtx.X() - tag_end.X(),2) + 
                             pow(mc_vtx.Y() - tag_end.Y(),2) + 
                             pow(mc_vtx.Z() - tag_end.Z(),2) );  

       if ( dist_st < 25 || dist_end < 25 ){
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
    else {
      //std::cout<<"\nEvent is : "<<_event <<", mult: "<<trk_map.size()<<std::endl ; 
      //If no true tracks aligned with reco track, mark it as noise
      _n_cosmic++;
      _bkgd_id = 2 ;
      _tree->Fill();

     return false;

    }

    auto mc_vtx = ev_mctrk->at(mc_max_it).Start() ;
    auto mc_end = ev_mctrk->at(mc_max_it).End() ;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Now have only events with true tagged tracks.  Next check containment of vertex
    ///////////////////////////////////////////////////////////////////////////////////////////////

    bool infv = true;

    if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
      infv = false;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// Now count the number of backgrounds and signals
    ///////////////////////////////////////////////////////////////////////////////////////////////
    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0;
    int n_mu = 0;
    int n_gamma = 0;

    // It turns out that for MC files, origin 2 is cosmic, 1 is neutrino. Can check this using 
    // BNB Only files or OpenCosmic files. The opposite is the case for reco tracks through back tracker, for some reason.
    if( ev_mctrk->at(mc_max_it).Origin() == 2 ){
      //std::cout<<"\nEvent is : "<<_event <<", "<<mc_max_dot<<std::endl ;
      _n_cosmic++;
      _bkgd_id = 2; 
    }
    else if( abs(nu.Nu().PdgCode()) == 12){
      _n_nue ++ ;
      _bkgd_id = 3;
    } 
    else if( nu.Nu().PdgCode() == -14){
      _n_antinumu++ ;
      _bkgd_id = 4;
    }
    else if( nu.Nu().PdgCode() == 14 && nu.CCNC() == 1 ){
      _n_nc++;
      _bkgd_id = 5;
    }

    if( _bkgd_id == -1 ){
      for ( auto const & p : parts ){
     
        if( p.StatusCode() == 1 && p.PdgCode() == 111 )
          n_pi0 ++;
        if( p.StatusCode() == 1 && p.PdgCode() == 13 )
          n_mu ++;
        if( p.StatusCode() == 1 && p.PdgCode() == 22)
          n_gamma ++;

      }   

        if( n_mu == 1){
          if( n_pi0 > 1 ){
            _bkgd_id = 6; 
            _multpi0++;
            //_n_ccother++;
          }
          else if( n_pi0 == 1 && !infv ){
            _bkgd_id = 7;
            _ccpi0_outfv ++ ;
	    //_n_ccother++;
          }
          else if( n_pi0 == 1 && infv ){ 
            _bkgd_id = 8;
            _tot_ccpi0 ++; 
            _event_list.emplace_back(_event);
          }
          else if( n_pi0 == 0 && n_gamma > 1 ){
            _bkgd_id = 9;
            _n_gammas++;
	    //_n_ccother++; 
          }
          else{
            _bkgd_id = 10;
            _n_ccother ++;
	    _ccother_list.emplace_back(_event);
          }
        }
      }
    
    _tree->Fill();    

    return true;
  }

  bool SingleGammaBackgroundCalc::finalize() {

    std::cout<<"Signals: "<<std::endl ;
    std::cout<<"Total CCpi0 : "<<_tot_ccpi0<<"/"<<_event_list.size()<<std::endl; 

    // Note that cc other includes secondary pi0s.
    std::cout<<"\nBackgrounds: "<<std::endl;
    std::cout<<"1) Noise : "<<_n_noise<< std::endl;
    std::cout<<"2) BNB Cosmic : "<<_n_cosmic<< std::endl;
    std::cout<<"3) Nue : "<<_n_nue<<std::endl;
    std::cout<<"4) Antinumus: "<<_n_antinumu<<std::endl;
    std::cout<<"5) NC pi0 : "<<_n_nc<<std::endl; // nc_pi0<<std::endl;
    std::cout<<"6) Multpi0s: "<<_multpi0<<std::endl ;
    std::cout<<"7) Nu outsude FV: "<<_ccpi0_outfv <<std::endl ;
    std::cout<<"8) Ngamma Events: "<<_n_gammas<<std::endl ;
    std::cout<<"9) OTHER CC Events: "<<_n_ccother<<std::endl ;

    std::cout<<"Total accounted backgrounds: "<< _n_cosmic + _n_nue + _n_antinumu + _n_nc + _multpi0 + _ccpi0_outfv + _n_gammas + _n_ccother <<std::endl ;
//    std::cout<<"Total accounted backgrounds: "<< _n_cosmic + _n_nue + _n_antinumu + _n_nc + _n_ccother <<std::endl ;

    //std::cout<<"\n\n"<<_event_list.size()<<" in Event list :" <<std::endl ;
    //for( auto const & e : _event_list) std::cout<<e<<", ";

    //std::cout<<"\n\n"<<_ccother_list.size()<<" in CCOther list :" <<std::endl ;
    //for( auto const & e : _ccother_list) std::cout<<e<<", ";

    if ( _fout ){
      _fout->cd();
      _tree->Write();
    }
  
    return true;
  }

}
#endif
