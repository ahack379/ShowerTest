#ifndef LARLITE_JUSTIFYPI0CUTS_CXX
#define LARLITE_JUSTIFYPI0CUTS_CXX

#include "JustifyPi0Cuts.h"
#include "DataFormat/shower.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/cluster.h"
#include "DataFormat/pfpart.h"

namespace larlite {

  bool JustifyPi0Cuts::initialize() {
    
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

   std::cout<<"\n\n****** "<<_event_list.size()<<" events found by 2 Shower Pi0Reco Module! ******"<<std::endl; 
   for ( auto const & e : _event_list )
     std::cout<<e <<", "; 

  std::cout<<"\n\n\n";
  
    return true;
  }

}
#endif
