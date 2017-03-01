#ifndef LARLITE_TRACKSTITCHER_CXX
#define LARLITE_TRACKSTITCHER_CXX

#include "TrackStitcher.h"
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

  bool TrackStitcher::initialize() {

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

     _tree->Branch("vtx_diff",&_vtx_diff,"vtx_diff/F"); 
     _tree->Branch("vtx_diff",&_vtx_diff,"vtx_diff/F"); 
    }   

    return true;
  }
  
  bool TrackStitcher::analyze(storage_manager* storage) {

    //auto const& geomH = ::larutil::GeometryHelper::GetME();
    // Track stitching plan: 
    // 0) Find events with multiplicity 1 
    // 1) Calculate the dot product between st,end of that track
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

    auto ev_trk = storage->get_data<event_track>("pandoraNu");
    if ( !ev_trk || !ev_trk->size() ) return false; 

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



    auto vtx = ev_vtx->at(0); 
    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };

    auto traj = nu.Trajectory().at(nu.Trajectory().size()-1);
    _nu_energy = traj.E(); 

    std::vector<int> trk_ids;
    trk_ids.reserve(ev_trk->size());
    
    //Map of lengths -> track id
    std::multimap<float,int> trk_map ;
    _length_v.clear();

    if( _ccpi0_event != 1 && _nue_event != 1) return false ;
    //std::cout<<"\nNew event! "<<_event<<std::endl ;

    int rad = 6;

    // Find closest + longest pandoraNu track to vertex  
    //std::cout<<"NUMBER TRACKS :" <<ev_trk->size()<<std::endl;
    int track_it = -1;
    int traj_size = 0;
    bool startClosest = true ;

    for ( size_t ti = 0; ti < ev_trk->size(); ti++ ) { 

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
       if ( dist_st < rad || dist_end < rad ){
          
	  if ( dist_st > rad ) startClosest = false ;

          _length_v.emplace_back(t_length);
          trk_ids.emplace_back(ti);
          trk_map.emplace(1./t_length,ti);

          // Quantities only used for single track mult study
          _length = t_length ; 
          track_it = ti ;
	  traj_size = ev_trk->at(ti).NumberTrajectoryPoints();
           }

        }

       //std::cout<<"Tracks within "<<rad<<" cm: "<<_length_v.size()<<std::endl; 
       _track_mult = _length_v.size() ;

      if(_track_mult == 1 && _ccpi0_event == 1) 
        std::cout<<"***************CCPi0 Event is: "<<_event<<", "<<_length<<", "<<std::endl;

      if(_track_mult == 1 && _nue_event == 1) 
        std::cout<<"****************Nue Event is: "<<_event<<", "<<_length<<", "<<std::endl;

      std::multimap<float,std::pair<int,int>> dot_to_trk_ids ;

      if( _track_mult == 1 && traj_size > 10 ){

            auto t_vtx_i = ev_trk->at(track_it).Vertex() ;
            auto t_end_i = ev_trk->at(track_it).End() ;
            auto dir = t_end_i - t_vtx_i ;
	    std::cout<<"Location of st: "<<t_vtx_i.X()<<", "<<t_vtx_i.Y()<<", "<<t_vtx_i.Z()<<std::endl ;
	    std::cout<<"Location of end: "<<t_end_i.X()<<", "<<t_end_i.Y()<<", "<<t_end_i.Z()<<std::endl ;

            auto t_end_use = t_end_i ; 
	    float t_length_i  = 0;
            TVector3 t_near_end, t_mom ;

            if ( !startClosest ){

	       t_end_use = t_vtx_i ;

               ev_trk->at(track_it).TrajectoryAtPoint(traj_size - 5 - 1,t_near_end,t_mom);
	       dir = t_near_end - t_end_i ;
               t_length_i = sqrt(pow(t_near_end.X() - t_end_i.X(),2) + pow(t_near_end.Y() - t_end_i.Y(),2) 
                                  +pow(t_near_end.Z() - t_end_i.Z(),2) );
	      }
            else {

               ev_trk->at(track_it).TrajectoryAtPoint(5,t_near_end,t_mom);
	       dir = t_near_end - t_vtx_i ;
               t_length_i = sqrt(pow(t_near_end.X() - t_vtx_i.X(),2) + pow(t_near_end.Y() - t_vtx_i.Y(),2) 
                                  +pow(t_near_end.Z() - t_vtx_i.Z(),2) );
	       }

            for ( size_t j = 0; j < ev_trk->size(); j++ ) { 

	      int traj_size2 = ev_trk->at(j).NumberTrajectoryPoints() ;
 
              if ( j == track_it || traj_size2 < 10 ) continue;
            
              auto t_vtx = ev_trk->at(j).Vertex() ;
              auto t_end = ev_trk->at(j).End() ;
              auto dir2 = t_vtx - t_end ;

              float dist_st = sqrt(pow(t_vtx.X() - t_end_use.X(),2) + pow(t_vtx.Y() - t_end_use.Y(),2) 
                                   +pow(t_vtx.Z() - t_end_use.Z(),2) );

              float dist_end = sqrt(pow(t_end.X() - t_end_use.X(),2) + pow(t_end.Y() - t_end_use.Y(),2) 
                                    +pow(t_end.Z() - t_end_use.Z(),2) );

              float dist_less = dist_st < dist_end ? dist_st : dist_end ;
              bool startClosest_j = dist_st < dist_end ? 1 : 0 ;

              float t_length = 0; 
	      TVector3 t_near ;

	      if ( startClosest_j ){

                ev_trk->at(j).TrajectoryAtPoint(5,t_near,t_mom);
	        dir2 = t_near - t_vtx ;
                t_length = sqrt(pow(t_near.X() - t_vtx.X(),2) + pow(t_near.Y() - t_vtx.Y(),2) 
                                  +pow(t_near_end.Z() - t_vtx.Z(),2) );
				  }
	      else{
                ev_trk->at(j).TrajectoryAtPoint(traj_size2 - 5 - 1,t_near,t_mom);
	        dir2 = t_near - t_end ;
                t_length = sqrt(pow(t_near.X() - t_end.X(),2) + pow(t_near.Y() - t_end.Y(),2) 
                                  +pow(t_near_end.Z() - t_end.Z(),2) );
	          }
	        
	      auto norm1 = sqrt( dir.X() * dir.X() + dir.Y() * dir.Y() + dir.Z() * dir.Z() );
	      auto norm2 = sqrt( dir2.X() * dir2.X() + dir2.Y() * dir2.Y() + dir2.Z() * dir2.Z() ) ;
	      //std::cout<<"Direction 2: "<< dir2.X()<<", "<<dir2.Y()<<", "<<dir2.Z()<<std::endl ;

              auto dot = 1. / norm1 / norm2 * (dir2.X() * dir.X() + dir2.Y() * dir.Y() + dir2.Z() * dir.Z() );
              //std::cout<<"Track entry "<<dot<<", "<<dist_less<<", "<<std::endl ;
	      //std::cout<<"TESTESTEST : "<<norm1 <<", "<<norm2<<", "<<dot<<", "<<dir.X()<<std::endl ;
              
              if (std::abs(dot) > 0.9 && dist_less < 8 ){
                //auto p = std::make_pair(i,j);
                //dot_to_trk_ids.emplace(dot,p);            
        	_length = t_length_i + t_length ;
                std::cout<<"Dot: "<<dot<<", "<<t_length_i<<", "<<t_length<<", "<<dist_less<<std::endl ;
               }
              } 
            }
    
    _tree->Fill();

    return true;
  }

  bool TrackStitcher::finalize() {

    if(_fout) { _fout->cd(); _tree->Write(); }

    return true;
  }

}
#endif
