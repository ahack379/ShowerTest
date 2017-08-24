#ifndef LARLITE_LOWLEVEL_CXX
#define LARLITE_LOWLEVEL_CXX

#include "LowLevel.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/hit.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/shower.h"

#include "LArUtil/GeometryHelper.h"

namespace larlite {

  bool LowLevel::initialize() {

    if(!_low_level_tree){
      _low_level_tree = new TTree("low_level_tree","");

      // Track info
      _low_level_tree->Branch("entry",&_entry,"entry/I");
      _low_level_tree->Branch("startx",&_startx,"startx/F");
      _low_level_tree->Branch("starty",&_starty,"starty/F");
      _low_level_tree->Branch("startz",&_startz,"startz/F");
      _low_level_tree->Branch("endx",&_endx,"endx/F");
      _low_level_tree->Branch("endy",&_endy,"endy/F");
      _low_level_tree->Branch("endz",&_endz,"endz/F");
      _low_level_tree->Branch("length",&_length,"length/F");
      _low_level_tree->Branch("theta",&_theta,"theta/F");
      _low_level_tree->Branch("phi",&_phi,"phi/F");

      // Vertex info
      _low_level_tree->Branch("vtxx",&_vtxx,"vtxx/F");
      _low_level_tree->Branch("vtxy",&_vtxy,"vtxy/F");
      _low_level_tree->Branch("vtxz",&_vtxz,"vtxz/F");
      _low_level_tree->Branch("vtxw",&_vtxw,"vtxw/F");
      _low_level_tree->Branch("vtxt",&_vtxt,"vtxt/F");
      _low_level_tree->Branch("vtx_trk_dist",&_vtx_trk_dist,"vtx_trk_dist/F");
      _low_level_tree->Branch("vtx_mc_reco_dist",&_vtx_mc_reco_dist,"vtx_mc_reco_dist/F");

      _low_level_tree->Branch("mc_vtxx",&_mc_vtxx,"mc_vtxx/F");
      _low_level_tree->Branch("mc_vtxy",&_mc_vtxy,"mc_vtxy/F");
      _low_level_tree->Branch("mc_vtxz",&_mc_vtxz,"mc_vtxz/F");
      _low_level_tree->Branch("mc_vtxw",&_mc_vtxw,"mc_vtxw/F");
      _low_level_tree->Branch("mc_vtxt",&_mc_vtxt,"mc_vtxt/F");
    
      // Flash info is not saved :(
    }

    if(!_hit_tree){
      _hit_tree = new TTree("hit_tree","");
      _hit_tree->Branch("entry",&_entry,"entry/I");
      _hit_tree->Branch("plane",&_plane,"plane/I");
      _hit_tree->Branch("charge",&_charge,"charge/F");
      _hit_tree->Branch("time_peak",&_time_peak,"time_peak/F");
      _hit_tree->Branch("time_width",&_time_width,"time_width/F");
      _hit_tree->Branch("gof",&_gof,"gof/F");
    }

    if(!_shower_tree){
      _shower_tree = new TTree("shr_tree","");
      _shower_tree->Branch("entry",&_entry,"entry/I");
      _shower_tree->Branch("shr_startx",&_shr_startx,"shr_startx/F");
      _shower_tree->Branch("shr_starty",&_shr_starty,"shr_starty/F");
      _shower_tree->Branch("shr_startz",&_shr_startz,"shr_startz/F");
      _shower_tree->Branch("shr_startw",&_shr_startw,"shr_startw/F");
      _shower_tree->Branch("shr_startt",&_shr_startt,"shr_startt/F");
      _shower_tree->Branch("shr_dirx",&_shr_dirx,"shr_dirx/F");
      _shower_tree->Branch("shr_diry",&_shr_diry,"shr_diry/F");
      _shower_tree->Branch("shr_dirz",&_shr_dirz,"shr_dirz/F");
      _shower_tree->Branch("shr_energy",&_shr_energy,"shr_energy/F");
      _shower_tree->Branch("shr_oangle",&_shr_oangle,"shr_oangle/F");
      _shower_tree->Branch("shr_dedx",&_shr_dedx,"shr_dedx/F");
      _shower_tree->Branch("shr_vtx_dist",&_shr_vtx_dist,"shr_vtx_dist/F");
      _shower_tree->Branch("shr_trk_delta_theta",&_shr_trk_delta_theta,"shr_trk_delta_theta/F");
      _shower_tree->Branch("shr_trk_delta_phi",&_shr_trk_delta_phi,"shr_trk_delta_phi/F");
    }

    _entry = -1;

    return true;
  }

  void LowLevel::clear(){
  
    _nhits0 = 0;
    _nhits1 = 0;
    _nhits2 = 0;

    _plane = -1;
    _charge = 0;
    _wire = -999;
    _time_peak = 0;
    _time_width = 0;
    _gof = -999;

    _startx = -999 ;
    _starty = -999 ;
    _startz = -999 ;
    _endx = -999;
    _endy = -999;
    _endz = -999;
    _length = -999;
    _theta = -999;
    _phi = -999;

    _vtxx = -999;
    _vtxy = -999;
    _vtxz = -999;
    _vtxw = -999;
    _vtxt = -999;
    _vtx_trk_dist = -999;
    _vtx_mc_reco_dist = -999;

    _mc_vtxx = -999;
    _mc_vtxy = -999;
    _mc_vtxz = -999;
    _mc_vtxw = -999;
    _mc_vtxt = -999;

    _shr_startx = -999;
    _shr_starty = -999;
    _shr_startz = -999;
    _shr_startw = -999;
    _shr_startt = -999;
    _shr_dirx = -999;
    _shr_diry = -999;
    _shr_dirz = -999;
    _shr_energy = -999;
    _shr_oangle = -999;
    _shr_dedx = -999;
    _shr_vtx_dist = -999;
    _shr_trk_delta_theta = -999;
    _shr_trk_delta_phi = -999;

  }
  
  bool LowLevel::analyze(storage_manager* storage) {

    clear();
    _entry++ ;

    auto geomH = ::larutil::GeometryHelper::GetME();

    // Store some truth info
    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) { 
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
    }   

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    auto traj = nu.Nu().Trajectory();
    _mc_vtxx = traj.at(traj.size() - 1).X();
    _mc_vtxy = traj.at(traj.size() - 1).Y();
    _mc_vtxz = traj.at(traj.size() - 1).Z();

    std::vector<float> mc_to_proj = { _mc_vtxx, _mc_vtxy, _mc_vtxz } ;
    auto mc2d = geomH->Point_3Dto2D(mc_to_proj,2) ;
    _mc_vtxw = mc2d.w ; 
    _mc_vtxt = mc2d.t ; 
 
    // Store hit info
    auto ev_hit = storage->get_data<event_hit>("gaushit");

    if ( ev_hit && ev_hit->size() ){
      
      for ( auto const & h : *ev_hit ){
        _plane = h.WireID().Plane ;
        _charge = h.Integral() ;
        _wire = h.WireID().Wire ;
        _time_peak = h.PeakTime();
        _time_width = h.SigmaPeakTime() ;
        _gof = h.GoodnessOfFit() ;

        _hit_tree->Fill(); 
      }
    }

    // Store track + vtx info
    auto ev_trk = storage->get_data<event_track>("numuCC_track");
    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");

    if ( ev_trk && ev_trk->size() && ev_vtx->size() ){
      auto t = ev_trk->at(0) ;
      auto v = ev_vtx->at(0) ;

      _startx = t.Vertex().X() ;
      _starty = t.Vertex().Y() ;
      _startz = t.Vertex().Z() ;
      _endx = t.End().X() ;
      _endy = t.End().Y() ;
      _endz = t.End().Z() ;
      _length = t.Length() ;
      //_length = sqrt( pow(startx - endx,2) + pow(starty - endy,2) + pow(startz - endz,2) ); 
      _theta = t.Theta() ;
      _phi = t.Phi() ; 

      std::vector<double> dir = { (_endx - _startx) / _length,
                                  (_endy - _starty) / _length,
                                  (_endz - _startz) / _length };

      auto dir_start = t.VertexDirection();
      std::vector<double> other_dir = { dir_start.X(), dir_start.Y(), dir_start.Z() }; 

      float dotProd = dir.at(0) * other_dir.at(0) + dir.at(1) * other_dir.at(1) +  dir.at(2) * other_dir.at(2) ;

      if( dotProd < 0 ) {
         TVector3 new_dir(-dir_start.X(),-dir_start.Y(),-dir_start.Z());
         _theta = new_dir.Theta();
         _phi = new_dir.Phi();
      }
      
      _vtxx = v.X();
      _vtxy = v.Y();
      _vtxz = v.Z();

      std::vector<float> reco_to_proj = { _vtxx, _vtxy, _vtxz } ;
      auto reco2d = geomH->Point_3Dto2D(reco_to_proj,2) ;
      _vtxw = reco2d.w ; 
      _vtxt = reco2d.t ; 
      
      auto st_dist = sqrt( pow(_vtxx - _startx,2) +pow(_vtxy - _starty,2) +pow(_vtxz - _startz,2)); 
      auto end_dist = sqrt( pow(_vtxx - _endx,2) +pow(_vtxy - _endy,2) +pow(_vtxz - _endz,2)); 

      _vtx_trk_dist = st_dist < end_dist ? st_dist : end_dist ;
      _vtx_mc_reco_dist = sqrt( pow(_mc_vtxx - _vtxx,2) + pow(_mc_vtxy - _vtxy,2) + pow(_mc_vtxz - _vtxz,2) );

     }

    _low_level_tree->Fill();

    auto ev_shr = storage->get_data<event_shower>("showerreco");

    if ( ev_shr && ev_shr->size() ){

      for( auto const & s : *ev_shr ){
        _shr_startx = s.ShowerStart().X();
        _shr_starty = s.ShowerStart().Y();
        _shr_startz = s.ShowerStart().Z();

        std::vector<float> shr_to_proj = { _shr_startx, _shr_starty, _shr_startz } ;
        auto shr2d = geomH->Point_3Dto2D(shr_to_proj,2) ;
        _shr_startw = shr2d.w ;
        _shr_startt = shr2d.t ;

        _shr_dirx = s.Direction().X();
        _shr_diry = s.Direction().Y();
        _shr_dirz = s.Direction().Z();

        _shr_dedx = s.Energy(2);
        _shr_oangle = s.OpeningAngle();
        _shr_dedx = s.dEdx(2);

        _shr_vtx_dist = sqrt( pow(_vtxx - _shr_startx,2) +
                              pow(_vtxy - _shr_starty,2) +
                              pow(_vtxz - _shr_startz,2) ); 
        
        _shr_trk_delta_theta = s.Direction().Theta() - _theta ;
        _shr_trk_delta_phi = s.Direction().Phi() - _phi ;

        _shower_tree->Write() ;
      }
    }

    return true;
  }

  bool LowLevel::finalize() {

    if ( _fout ){
      _fout->cd();
      _low_level_tree->Write();
      _hit_tree->Write();
      _shower_tree->Write() ;
    }
     
    return true;
  }

}
#endif
