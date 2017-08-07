#ifndef LARecoRLITE_RECORLENERGYCUTSTUDY_CXX
#define LARecoRLITE_RECORLENERGYCUTSTUDY_CXX

#include "RecoRLEnergyCutStudy.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"

namespace larlite {

  bool RecoRLEnergyCutStudy::initialize() {

    if( !_shower_tree ){
     _shower_tree = new TTree("shower_tree","");
     _shower_tree->Branch("event",&_event,"event/I");
     _shower_tree->Branch("nshrs",&_nshrs,"nshrs/I");
     _shower_tree->Branch("signal",&_signal,"signal/B");
     _shower_tree->Branch("shower_e",&_shower_e,"shower_e/F");
     _shower_tree->Branch("shower_rl",&_shower_rl,"shower_rl/F");
     _shower_tree->Branch("shower_angle",&_shower_angle,"shower_angle/F");
     _shower_tree->Branch("shower_oangle",&_shower_oangle,"shower_oangle/F");
     _shower_tree->Branch("mu_angle",&_mu_angle,"mu_angle/F");
    }

    if( !_hit_tree){
     _hit_tree= new TTree("hit_tree","");
     _hit_tree->Branch("event",&_event,"event/I");
     _hit_tree->Branch("nshrs",&_nshrs,"nshrs/I");
     _hit_tree->Branch("signal",&_signal,"signal/B");
     _hit_tree->Branch("shower_e",&_shower_e,"shower_e/F");
     _hit_tree->Branch("shower_rl",&_shower_rl,"shower_rl/F");
     _hit_tree->Branch("shower_angle",&_shower_angle,"shower_angle/F");
     _hit_tree->Branch("shower_oangle",&_shower_oangle,"shower_oangle/F");
     _hit_tree->Branch("mu_angle",&_mu_angle,"mu_angle/F");
     _hit_tree->Branch("hit_ratio",&_hit_ratio,"std::vector<float>");
     _hit_tree->Branch("gaus_hits_per_r","std::vector<float>",&_gaus_hits_per_r);
     _hit_tree->Branch("shr_hits_per_r","std::vector<float>",&_shr_hits_per_r);
     _hit_tree->Branch("radii","std::vector<float>",&_radii);
    }

    _event = -1;
    _signal = false; 

    return true;
  }

 void RecoRLEnergyCutStudy::clear(){

    _shower_rl = -10;
    _shower_e  = -10;
    _shower_angle  = -10;
    _shower_oangle  = -10;
    _mu_angle  = -10;
    //_hit_ratio = -10 ;

    _signal = false ;

    _nshrs = 0;

  }

  bool RecoRLEnergyCutStudy::analyze(storage_manager* storage) {

    _event++; 
    clear();

    int rad_its = 10;

    _radii.clear();
    _hit_ratio.clear();
    _gaus_hits_per_r.clear();
    _shr_hits_per_r.clear();

    _radii.reserve(rad_its);
    _hit_ratio.reserve(rad_its);
    _gaus_hits_per_r.reserve(rad_its);
    _shr_hits_per_r.reserve(rad_its);

    std::cout<<"\nEvent : "<<_event <<std::endl;

    auto ev_mctruth = storage->get_data<event_mctruth>("generator");
    if( !ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Not enough mctruth..." <<ev_mctruth->size()<<std::endl; return false; }

    auto ev_s = storage->get_data<event_shower>("showerreco");
    if( !ev_s || !ev_s->size() ) {
      std::cout<<"Not enough showers..." <<ev_s->size()<<std::endl; return false; }

    //auto ev_t = storage->get_data<event_track>("numuCC_track");
    //if( !ev_t || !ev_t->size() ) {
    //  std::cout<<"No track???..." <<ev_t->size()<<std::endl; return false; }

    //auto tag_trk = ev_t->at(0) ;
    //_mu_angle = cos(tag_trk.Theta());

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    if( xyz[0] < 0 || xyz[0] > 256.35 || xyz[1] > 116.5 || xyz[1] < -116.5 || xyz[2] < 0 || xyz[2] > 1036.8 )
      return false;

    auto ev_hit_g = storage->get_data<event_hit>("gaushit");
    if ( !ev_hit_g || !ev_hit_g->size() ) {std::cout<<"Returning..."<<std::endl ; return false; }

    //for(int j = 0; j < rad_its; j++){

    //  int rad = float ( (j+1) * 5. ) ;
    //  int hits_in_rad = 0.;
    //  int hits_in_rad_g = 0.;

    //  for(auto const & h : *ev_hit_g){

    //    if( h.WireID().Plane != 2 ) continue;

    //    auto w = h.WireID().Wire * geomH->WireToCm();
    //    auto t = h.PeakTime() * geomH->TimeToCm() ;

    //    auto dist = sqrt( pow(w - vtx_w,2) + pow(t - vtx_t,2) );

    //    if(dist <= rad) hits_in_rad_g ++ ;
    //    if(dist <= rad && h.GoodnessOfFit() >= 0) hits_in_rad++ ;

    // }

    // _radii.emplace_back(rad);
    // _shr_hits_per_r.emplace_back(hits_in_rad);
    // _gaus_hits_per_r.emplace_back(hits_in_rad_g);

    // if( hits_in_rad_g < 15 )
    //   _hit_ratio.emplace_back(0.);
    // else
    //   _hit_ratio.emplace_back(float(hits_in_rad)/hits_in_rad_g);
    //}

    auto parts = truth.GetParticles() ;

    int n_pi0 = 0;
    int n_mu = 0;

    for ( auto const & p : parts ){
      if ( p.PdgCode() == 111 && p.StatusCode() == 1 )
        n_pi0++;

      if ( p.PdgCode() == 13 && p.StatusCode() == 1 )
        n_mu++;
    }

    if ( n_pi0 == 1 && n_mu == 1 ) _signal = true;
  
    _nshrs = ev_s->size() ;

    for ( int i = 0; i < ev_s->size(); i++ ){ 

      auto s = ev_s->at(i); 

      ::geoalgo::Point_t vertex(xyz[0],xyz[1],xyz[2]);

      ::geoalgo::Vector_t rev(-1.*s.Direction()) ;
      auto s_bkwrd_hl = ::geoalgo::HalfLine_t(s.ShowerStart(),rev);

      auto IP = _geoAlgo.SqDist(vertex, s_bkwrd_hl) ;
      if( IP > 4 ) continue;

      auto st = s.ShowerStart() ;
      auto dist = sqrt( pow(st.X() - xyz[0],2) + pow(st.Y() - xyz[1],2) + pow(st.Z() - xyz[2],2) ); 
    
      _shower_e = s.Energy(2) ;
      _shower_rl = dist;
      _shower_angle = s.Direction().Y();
      _shower_oangle = s.OpeningAngle() ;
      

      _shower_tree->Fill();
    }  

    return true;
  }

  bool RecoRLEnergyCutStudy::finalize() {

    if(_fout) {
      _fout->cd();
      _shower_tree->Write();
    }
  
    return true;
  }

}
#endif
