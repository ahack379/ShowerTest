#ifndef LARLITE_NOMUONVTXDENSITY_CXX
#define LARLITE_NOMUONVTXDENSITY_CXX

#include "NoMuonVtxDensity.h"
#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include "LArUtil/GeometryHelper.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>

namespace larlite {

  bool NoMuonVtxDensity::initialize() {

   _event = 0 ;
   _hits_in_rad_g = 0;
   _hits_in_rad = 0;
   _hits_in_rad_ass = 0;

   if( !_tree ){
     _tree = new TTree("tree","tree");
     _tree->Branch("hit_ratio_v","std::vector<float>",&_hit_ratio_v); 
    }   

    return true;
  }
  
  bool NoMuonVtxDensity::analyze(storage_manager* storage) {


    auto const& geomH = ::larutil::GeometryHelper::GetME();

    std::cout<<"\nNew event! "<<_event<<std::endl ;
    _event++;

    int rad_its = 15;
    _hit_ratio_v.clear();
    _hit_ratio_v.reserve(rad_its);
    
    auto ev_hit = storage->get_data<event_hit>("hit02");
    if ( !ev_hit || !ev_hit->size() ) return false;

    auto ev_hit_g = storage->get_data<event_hit>("gaushit");
    if ( !ev_hit_g || !ev_hit_g->size() ) return false; 

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");
    if ( !ev_vtx || !ev_vtx->size() ) return false; 

    auto ev_trk = storage->get_data<event_track>("pandoraNu");
    if ( !ev_trk || !ev_trk->size() ) return false; 

    auto vtx = ev_vtx->at(0); 
    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };

    std::vector<int> trk_ids;
    trk_ids.reserve(ev_trk->size());
    
    //Map of lengths -> track id
    std::multimap<float,int> trk_map ;

    // Find closest + longest pandoraNu track to vertex
    for ( size_t ti = 0; ti < ev_trk->size(); ti++ ) { //auto const & v : ev_vtx_pandoraNu ){

      auto t_vtx = ev_trk->at(ti).Vertex() ;
      auto t_end = ev_trk->at(ti).End() ;
       
      float dist_st = sqrt( pow(t_vtx.X() - vtxXYZ[0],2) + 
                            pow(t_vtx.Y() - vtxXYZ[1],2) + 
                            pow(t_vtx.Z() - vtxXYZ[2],2) ); 

      float dist_end = sqrt( pow(t_end.X() - vtxXYZ[0],2) + 
                             pow(t_end.Y() - vtxXYZ[1],2) + 
                             pow(t_end.Z() - vtxXYZ[2],2) ); 
       if ( dist_st < 2 || dist_end < 2 ){
          float len = ev_trk->at(ti).Length();
          trk_ids.emplace_back(ti);
          trk_map.emplace(1./len,ti);
           }
      }
    //std::cout<<"Longest track + track ID : "<<1./trk_map.begin()->first<<", "<<trk_map.begin()->second<<std::endl ;
    //std::cout<<trk_ids.size()<<" tracks near vertex "<<std::endl ;

    // Want to ignore longest track associated to vertex-- Get handle to association 
    auto ev_ass = storage->get_data<larlite::event_ass>("pandoraNu");

    auto ev_hit_cosRem = storage->get_data<event_hit>("pandoraCosmicKHitRemoval");

    // Get association to trk => hit and hit => trk
    auto const& ass_hit_v = ev_ass->association(ev_trk->id(), ev_hit_cosRem->id());

    if (!ass_hit_v.size()) {
      std::cout << "No ass from track => hit! " << std::endl;
      return false;
      }
    
    //std::cout<<"Gaus, Shr, Ass hits: "<<ev_hit_g->size()<<", "<<ev_hit->size()<<", "<<ass_hit_v.at(trk_map.begin()->second).size()<<std::endl ;

    //std::cout<<"Hit ass size at track 1! "<<ass_hit_v.at(1).size() <<std::endl ;

    auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,2);
    
    auto vtx_w = vtxWT.w; 
    auto vtx_t = vtxWT.t + 800 * geomH->TimeToCm() ; 

    float rad = 0. ;

    std::vector<float> x(50,0);
    std::vector<float> y(50,0);

    for(int j = 0; j < rad_its; j++){

      rad = float ( (j+1) * 5. ) ; 
      _hits_in_rad = 0.; 
      _hits_in_rad_g = 0.; 
      _hits_in_rad_ass = 0;
      std::vector<::cv::Point> contour;

      // Build circle around radius
      for(int i = 0; i < x.size(); i++){

         x[i] = vtx_w + rad * cos(M_PI * 2. * i / x.size());
         y[i] = vtx_t + rad * sin(M_PI * 2. * i / x.size());

         ::cv::Point pt(x[i],y[i]) ;
         contour.emplace_back(pt);
         }

       ::cv::Point pt(contour.at(0).x, contour.at(0).y) ;
       contour.emplace_back(pt);

      // Count shower hits in radius
      for(int i = 0; i < ev_hit_g->size(); i++){

          auto const & h_g = ev_hit_g->at(i); 

          if( h_g.WireID().Plane == 2 ){ 
       
            ::cv::Point h_pt_g(h_g.WireID().Wire * geomH->WireToCm(),
                             h_g.PeakTime() * geomH->TimeToCm()) ;
            auto inside_g = ::cv::pointPolygonTest(contour,h_pt_g,false);

            if(inside_g >= 0) _hits_in_rad_g++ ;
            }
          
          if( i < ev_hit->size() ){

              auto const & h = ev_hit->at(i);
              
              if( h.WireID().Plane == 2 ){
                
                ::cv::Point h_pt(h.WireID().Wire * geomH->WireToCm(),
                                 h.PeakTime() * geomH->TimeToCm()) ;
                auto inside = ::cv::pointPolygonTest(contour,h_pt,false);

                if(inside  >= 0) _hits_in_rad ++ ;
                }
              }


           if( trk_map.size() ){ 
             if( i < ass_hit_v.at(trk_map.begin()->second).size() ){ 

                auto h_i = ass_hit_v.at(trk_map.begin()->second).at(i) ;
                auto const & h_ass = ev_hit_cosRem->at(h_i);

                 if( h_ass.WireID().Plane != 2 ) continue;

                ::cv::Point h_pt_ass(h_ass.WireID().Wire * geomH->WireToCm(),
                                 h_ass.PeakTime() * geomH->TimeToCm()) ;

                auto inside_ass = ::cv::pointPolygonTest(contour,h_pt_ass,false);

                if(inside_ass >= 0) _hits_in_rad_ass++ ;
              }
	    }
        }

       //std::cout<<"Gauss hits, Ass hits, Shr hits: "<<_hits_in_rad_g<<", "<<_hits_in_rad_ass<<", "<<_hits_in_rad<<std::endl ;

      _hits_in_rad_g -= _hits_in_rad_ass ;
      //_hits_in_rad -= _hits_in_rad_ass ;

      if( _hits_in_rad_g <= 0 ) 
        _hit_ratio_v.emplace_back(0.);
      else 
        _hit_ratio_v.emplace_back(float(_hits_in_rad)/_hits_in_rad_g);
      //std::cout<<"Hit ratio : "<<_hit_ratio_v.at(_hit_ratio_v.size() - 1)<<std::endl;
     }
  
    _tree->Fill();

    // if( _hits_in_rad_g == 0 ){
    //   //_hit_ratio.emplace_back(0.);
    //   return false ;
    //    }
    // else{
    //   float ratio = float( _hits_in_rad ) / _hits_in_rad_g ;
    //   //_hit_ratio.emplace_back( ratio );
    //    
    //   if( ratio < _ratio_cut ){
    //      return false; 
    //     }
    //   else 
    //      if( IsSignal(intType) )
    //        _good_tags ++ ;
    //      else
    //        _bad_tags ++ ;
    //      return true; 
    //    }

     //std::cout<<"Hits in rad "<<rad<<": "<<float(_hits_in_rad)/_hits_in_rad_g<<std::endl; 

    return true;
  }

  bool NoMuonVtxDensity::finalize() {

    if(_fout) { _fout->cd(); _tree->Write(); }

    //std::cout<<float(_good_tags)/_total_sig_events*100<<"\% of signal events in Selection 2 pass this filter"<<std::endl ;
    //std::cout<<float(_bad_tags)/(_good_tags + _bad_tags)*100<<"\% of our final sample is shit we don't want"<<std::endl ;
  
    return true;
  }


// t is interaction type
bool NoMuonVtxDensity::IsSignal(int t){

    if( (t == 1004 || t == 1011 || t == 1080 || t == 1086 || t == 1090) )  
      return true ;
    else 
      return false ;
  }

}
#endif
