#ifndef LARLITE_HITDENSITY_CXX
#define LARLITE_HITDENSITY_CXX

#include "HitDistDensity.h"
#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"
#include "LArUtil/GeometryHelper.h"

namespace larlite {

  bool HitDistDensity::initialize() {

   if( !_tree ){
     _tree = new TTree("tree","tree"); 
     _tree->Branch("hits_tot",&_hits_tot,"hits_tot/I"); 
     _tree->Branch("radii","std::vector<float>",&_radii); 
     _tree->Branch("event",&_event,"event/I"); 
     _tree->Branch("hits_per_r","std::vector<float>",&_hits_per_r); 
     _tree->Branch("gaus_hits_per_r","std::vector<float>",&_gaus_hits_per_r); 
     _tree->Branch("shr_hits_per_r","std::vector<float>",&_shr_hits_per_r); 

     _tree->Branch("shr_hits_tot",&_shr_hits_tot,"shr_hits_tot/I"); 

     _tree->Branch("vtx_x",&_vtx_x,"vtx_x/F"); 
     _tree->Branch("vtx_y",&_vtx_y,"vtx_y/F"); 
     _tree->Branch("vtx_z",&_vtx_z,"vtx_z/F"); 
    }

    _event = 0;
    _keep = 0;

    return true;
  }
  
  bool HitDistDensity::analyze(storage_manager* storage) {

    int rad_its = 15;
    _hits_tot = 0.;
    _shr_hits_tot = 0.;

    _radii.clear();
    _density.clear();
    _hits_per_r.clear();
    _gaus_hits_per_r.clear();
    _shr_hits_per_r.clear();

    _radii.reserve(rad_its);
    _density.reserve(rad_its);
    _hits_per_r.reserve(rad_its);
    _gaus_hits_per_r.reserve(rad_its);
    _shr_hits_per_r.reserve(rad_its);

    auto const& geomH = ::larutil::GeometryHelper::GetME();

    //std::cout<<"\nNew event! "<<_event<<", sub, run, event: "<<storage->subrun()<<", "<<storage->run()<<", "<<storage->event()<<std::endl ;
    _event++;

    auto ev_hit_g = storage->get_data<event_hit>("gaushit");//hit02");
    if ( !ev_hit_g || !ev_hit_g->size() ) {std::cout<<"Returning..."<<std::endl ; return false; }

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");
    if ( !ev_vtx || !ev_vtx->size() ) {std::cout<<"No vertex. Returning..."<<std::endl ; return false; }

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) return false;

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    // Check FV
    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();
    auto e = traj.at(traj.size() - 1).E();

    bool infv = true; 

    if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
        infv = false ;
 
    // Check that there is only 1 muon and only 1 pi0 originating from vertex
    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0; 
    int n_mu = 0;

    for ( auto const & p : parts ){

        if( p.StatusCode() == 1 && p.PdgCode() == 111 )
          n_pi0 ++ ;

        if( p.StatusCode() == 1 && p.PdgCode() == 13) 
          n_mu++ ;
    }

    if((n_mu != 1 || n_pi0 != 1 || !infv ) && _get_pi0s )
       return false;
    
    if((n_mu == 1 && n_pi0 == 1 && infv) && !_get_pi0s)
       return false;

    for( auto const & h : *ev_hit_g ){
      if ( h.WireID().Plane == 2 ) _hits_tot++;
      if ( h.WireID().Plane == 2 && h.GoodnessOfFit() >= 0 ) _shr_hits_tot++;
    }

    auto vtx = ev_vtx->at(0); 
    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };

    _vtx_x = vtx.X();
    _vtx_y = vtx.Y();
    _vtx_z = vtx.Z();

    auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,2);
    auto vtx_w = vtxWT.w ; // Comes in cm; 
    auto vtx_t = vtxWT.t + 800. * geomH->TimeToCm() ; // 800 for single particle files

    float rad = 0. ;

    std::vector<float> x(50,0);
    std::vector<float> y(50,0);

    for(int j = 0; j < rad_its; j++){

        rad = float ( (j+1) * 5. ) ;
        _hits_in_rad = 0.;
        _hits_in_rad_g = 0.;

        for(auto const & h : *ev_hit_g){

          if( h.WireID().Plane != 2 ) continue; 

          auto w = h.WireID().Wire * geomH->WireToCm();
          auto t = h.PeakTime() * geomH->TimeToCm() ;

          auto dist = sqrt( pow(w - vtx_w,2) + pow(t - vtx_t,2) );

          if(dist <= rad) _hits_in_rad_g ++ ;

          if(dist <= rad && h.GoodnessOfFit() >= 0) _hits_in_rad++ ;

	}

        _radii.emplace_back(rad);
	_density.emplace_back(_hits_in_rad / (M_PI * rad * rad )) ;
	_gaus_hits_per_r.emplace_back(_hits_in_rad_g);
	_shr_hits_per_r.emplace_back(_hits_in_rad);

        if( _hits_in_rad_g < 15 ) 
	  _hits_per_r.emplace_back(0.);
        else 
	  _hits_per_r.emplace_back(float(_hits_in_rad)/_hits_in_rad_g);
      }

    _tree->Fill();
  
    return true;
  }

  bool HitDistDensity::finalize() {

    //std::cout<<"At 60cm, for 0.24 cut we keep: "<<_keep<<std::endl ;

    if(_fout){
      _fout->cd(); 
      _tree->Write(); 
      }
  
    return true;
  }

}
#endif
