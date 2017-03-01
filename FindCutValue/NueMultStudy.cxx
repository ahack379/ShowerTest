#ifndef LARLITE_NUEMULTSTUDY_CXX
#define LARLITE_NUEMULTSTUDY_CXX

#include "NueMultStudy.h"
#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include "DataFormat/mctrack.h"
#include "LArUtil/GeometryHelper.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>

namespace larlite {

  bool NueMultStudy::initialize() {

   _event = -1 ;
   _track_mult = 0;

   if( !_tree ){
     _tree = new TTree("tree","tree");
     _tree->Branch("length_v","std::vector<float>",&_length_v); 
     _tree->Branch("track_mult",&_track_mult,"track_mult/I"); 
     _tree->Branch("nue_event",&_nue_event,"nue_event/I"); 
     _tree->Branch("ccpi0_event",&_ccpi0_event,"ccpi0_event/I"); 
     _tree->Branch("nu_energy",&_nu_energy,"nu_energy/F"); 
     _tree->Branch("length",&_length,"length/F"); 
     _tree->Branch("pid",&_pid,"pid/I"); 
     _tree->Branch("vtx_diff",&_vtx_diff,"vtx_diff/F"); 
    }   

    return true;
  }
  
  bool NueMultStudy::analyze(storage_manager* storage) {

    //auto const& geomH = ::larutil::GeometryHelper::GetME();
    // Track stitching plan: 
    // 0) Find events with multiplicity 1 
    // 1) Calculate the dot product between stst,stend,endst,endend of that track
    //    with other tracks in the event; do this dot product near the end of the track whose
    //    completeness is in question
    // 2) Look at dot product of these segments 
    

    _event++;
    _track_mult = 0;
    _nue_event = 0;
    _ccpi0_event = 0;
    _nu_energy = 0;
    _length = 0;
    _pid = 0;
    _vtx_diff = -1;

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");
    if ( !ev_vtx || !ev_vtx->size() ) return false; 

    //auto ev_mcvtx = storage->get_data<event_vertex>("mcvertex");
    //if ( !ev_mcvtx || !ev_mcvtx->size() ) return false; 

    auto ev_trk = storage->get_data<event_track>("pandoraNu");
    if ( !ev_trk || !ev_trk->size() ) return false; 

    //auto ev_trk = storage->get_data<event_mctrack>("mcreco");
    //if ( !ev_trk || !ev_trk->size() ) return false; 

    auto ev_truth = storage->get_data<event_mctruth>("generator");
    if ( !ev_truth || !ev_truth->size() ) return false; 
    
    auto const & nu = ev_truth->at(0).GetNeutrino().Nu();
    auto parts = ev_truth->at(0).GetParticles();
    int n_pi0 = 0;
    int n_mu = 0;

    for ( auto const & p : parts ){
    
      if( p.StatusCode() == 1 && p.PdgCode() == 111 )
        n_pi0 += 1;

      if( p.StatusCode() == 1 && p.PdgCode() == 13 )
        n_mu += 1;

      }   

    if( n_mu == 1 && n_pi0 == 1)
      _ccpi0_event = 1;

    if ( abs(nu.PdgCode()) == 12 )
      _nue_event = 1;


    auto vtx2 = nu.Trajectory().at(nu.Trajectory().size()-1);
    _nu_energy = vtx2.E(); 

    auto vtx = ev_vtx->at(0); 
    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };

    //auto mcvtx = ev_mcvtx->at(0); 
    //std::vector<double> mcvtxXYZ = { mcvtx.X(), mcvtx.Y(), mcvtx.Z() };
    //std::vector<double> vtxXYZ = { mcvtx.X(), mcvtx.Y(), mcvtx.Z() };

    //_vtx_diff = sqrt(pow(mcvtxXYZ[0] - vtxXYZ[0],2) +pow(mcvtxXYZ[1] - vtxXYZ[1],2) +pow(mcvtxXYZ[2] - vtxXYZ[2],2) );

    //if( vtxXYZ[0] < 20. || vtxXYZ[0] > 236.35 || vtxXYZ[1] < -96.5 || vtxXYZ[1] > 96.5 
    //|| vtxXYZ[2] < 10. || vtxXYZ[2] > 1026.8 )
     //   return false;

    std::vector<int> trk_ids;
    trk_ids.reserve(ev_trk->size());
    
    //Map of lengths -> track id
    std::multimap<float,int> trk_map ;
    _length_v.clear();

    if( _ccpi0_event != 1 && _nue_event != 1) return false ;
    //std::cout<<"\nNew event! "<<_event<<std::endl ;

    int j = 6;

      // Find closest + longest pandoraNu track to vertex  
      //std::cout<<"NUMBER TRACKS :" <<ev_trk->size()<<std::endl;
      int track_it = -1;

      for ( size_t ti = 0; ti < ev_trk->size(); ti++ ) { //auto const & v : ev_vtx_pandoraNu ){

        auto t_vtx = ev_trk->at(ti).Vertex(); //Start(); //Vertex() ;
        auto t_end = ev_trk->at(ti).End() ;
	auto t_length = sqrt(pow(t_vtx.X() - t_end.X(),2) + pow(t_vtx.Y() - t_end.Y(),2) 
	                    +pow(t_vtx.Z() - t_end.Z(),2) );
        
        float dist_st = sqrt( pow(t_vtx.X() - vtxXYZ[0],2) + 
                              pow(t_vtx.Y() - vtxXYZ[1],2) + 
                              pow(t_vtx.Z() - vtxXYZ[2],2) ); 

        float dist_end = sqrt( pow(t_end.X() - vtxXYZ[0],2) + 
                               pow(t_end.Y() - vtxXYZ[1],2) + 
                               pow(t_end.Z() - vtxXYZ[2],2) ); 

         //std::cout<<"Dist out: "<<dist_st<<", "<<dist_end <<std::endl;
         if ( dist_st < j || dist_end < j ){
            float len = t_length ; //ev_trk->at(ti).Length();
            _length_v.emplace_back(len);
            trk_ids.emplace_back(ti);
            trk_map.emplace(1./len,ti);

            // Quantities only used for single track mult study
	    _length = len ; 
	    track_it = ti ;
	    //_pid = ev_trk->at(ti).PdgCode() ;
	    //std::cout<<"length: "<<len<<std::endl ;
             }

          }

       //std::cout<<"Tracks within "<<j<<" cm: "<<_length_v.size()<<std::endl; 
       _track_mult = _length_v.size() ;

      if(_track_mult == 1 && _ccpi0_event == 1) 
        std::cout<<"***************CCPi0 Event is: "<<_event<<", "<<_vtx_diff<<std::endl;

      if(_track_mult == 1 && _nue_event == 1) 
        std::cout<<"****************Nue Event is: "<<_event<<", "<<_vtx_diff<<std::endl;

      std::multimap<float,std::pair<int,int>> dot_to_trk_ids ;

    _tree->Fill();

    return true;
  }

  bool NueMultStudy::finalize() {

    if(_fout) { _fout->cd(); _tree->Write(); }

    return true;
  }

}
#endif
