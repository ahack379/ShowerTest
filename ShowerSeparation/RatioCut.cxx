#ifndef LARLITE_RATIOCUT_CXX
#define LARLITE_RATIOCUT_CXX

#include "RatioCut.h"
#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"
#include "LArUtil/GeometryHelper.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>

namespace larlite {

  bool RatioCut::initialize() {

   if( !_tree ){
     _tree = new TTree("tree","tree"); 
     _tree->Branch("hits_in_rad",&_hits_in_rad,"hits_in_rad/F"); 
     _tree->Branch("hits_in_rad_g",&_hits_in_rad_g,"hits_in_rad_g/F"); 
     _tree->Branch("hits_tot",&_hits_tot,"hits_tot/F"); 
     _tree->Branch("event",&_event,"event/I"); 
     _tree->Branch("hits_per_r","std::vector<float>",&_hits_per_r); 
    }

    _event = 0;

    return true;
  }
  
  bool RatioCut::analyze(storage_manager* storage) {

    int rad_its = 15;
    _hits_tot = 0.;
//temp
    _hits_per_r.clear();
    _hits_per_r.reserve(rad_its);

    auto const& geomH = ::larutil::GeometryHelper::GetME();

    //std::cout<<"\nNew event! "<<_event<<std::endl ;
    _event++;

    auto ev_hit = storage->get_data<event_hit>("hit02"); //shrhits");//hit02");
    if ( !ev_hit || !ev_hit->size() ){std::cout<<"Returning..."<<std::endl ; return false; }

    auto ev_hit_g = storage->get_data<event_hit>("gaushit");//hit02");
    if ( !ev_hit_g || !ev_hit_g->size() ) {std::cout<<"Returning..."<<std::endl ; return false; }

    auto ev_vtx = storage->get_data<event_vertex>("mcvertex"); //numuCC_vertex");
    if ( !ev_vtx || !ev_vtx->size() ) {std::cout<<"Returning..."<<std::endl ; return false; }

    auto vtx = ev_vtx->at(0); 

    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };
    auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,2);//plane);
         
    auto vtx_w = vtxWT.w ; // Comes in cm; 
    auto vtx_t = vtxWT.t + 800. * geomH->TimeToCm() ; // 800 for single particle files

    std::vector<::cv::Point> contour;

    float rad = 60. ;
    float offset = 0.;

    std::vector<float> x(50,0);
    std::vector<float> y(50,0);

    _hits_in_rad = 0.;
    _hits_in_rad_g = 0.;

    // Build circle around radius
    for(int i = 0; i < x.size(); i++){

       x[i] = vtx_w + rad * cos(M_PI * 2. * i / x.size());
       y[i] = vtx_t + rad * sin(M_PI * 2. * i / x.size());

       ::cv::Point pt(x[i],y[i]) ;
       contour.emplace_back(pt);
       } 

     ::cv::Point pt(contour.at(0).x, contour.at(0).y) ;
     contour.emplace_back(pt);

     float ratio = 0.;
	

     // Count hits in radius
     for(auto const & h : *ev_hit){

         if( h.WireID().Plane != 2 ) continue; 

         ::cv::Point h_pt(h.WireID().Wire * geomH->WireToCm(),
                          h.PeakTime() * geomH->TimeToCm()) ;
         auto inside = ::cv::pointPolygonTest(contour,h_pt,false); 

         if(inside  >= 0) _hits_in_rad ++ ;
         }

     for(auto const & h : *ev_hit_g){

         if( h.WireID().Plane != 2 ) continue; 

         ::cv::Point h_pt(h.WireID().Wire * geomH->WireToCm(),
                          h.PeakTime() * geomH->TimeToCm()) ;
         auto inside = ::cv::pointPolygonTest(contour,h_pt,false); 

         if(inside  >= 0) _hits_in_rad_g ++ ;
         }

	//std::cout<<"At rad "<<rad<<" gaus and hit02 ratios: "<<_hits_in_rad_g <<", "<<_hits_in_rad <<std::endl ;

    if( ratio < 0.1 ) return false ;
    
    _hits_tot = ev_hit->size() ;
    _tree->Fill();
  
    return true;
  }

  bool RatioCut::finalize() {

    if(_fout){
      _fout->cd(); 
      _tree->Write(); 
      }
  
    return true;
  }

}
#endif
