#ifndef LARLITE_HITDENSITY_CXX
#define LARLITE_HITDENSITY_CXX

#include "VtxDensity.h"
#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"
#include "LArUtil/GeometryHelper.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>

namespace larlite {

  bool VtxDensity::initialize() {

   if( !_tree ){
     _tree = new TTree("tree","tree"); 
     //_tree->Branch("hits_in_rad",&_hits_in_rad,"hits_in_rad/F"); 
     //_tree->Branch("hits_in_rad_g",&_hits_in_rad_g,"hits_in_rad_g/F"); 
     _tree->Branch("hits_tot",&_hits_tot,"hits_tot/F"); 
     _tree->Branch("radii","std::vector<float>",&_radii); 
     //_tree->Branch("density","std::vector<float>",&_density); 
     _tree->Branch("event",&_event,"event/I"); 
     _tree->Branch("hits_per_r","std::vector<float>",&_hits_per_r); 
     _tree->Branch("gaus_hits_per_r","std::vector<float>",&_gaus_hits_per_r); 
     _tree->Branch("shr_hits_per_r","std::vector<float>",&_shr_hits_per_r); 
    }

    _event = 0;
    _keep = 0;

    return true;
  }
  
  bool VtxDensity::analyze(storage_manager* storage) {

    int rad_its = 15;
    _hits_tot = 0.;

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

    //auto ev_hit = storage->get_data<event_hit>("hit02"); //shrhits");//hit02");
    //if ( !ev_hit || !ev_hit->size() ){std::cout<<"Returning..."<<std::endl ; return false; }

    auto ev_hit_g = storage->get_data<event_hit>("gaushit");//hit02");
    if ( !ev_hit_g || !ev_hit_g->size() ) {std::cout<<"Returning..."<<std::endl ; return false; }

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");
    if ( !ev_vtx || !ev_vtx->size() ) {std::cout<<"No vertex. Returning..."<<std::endl ; return false; }

    _hits_tot = ev_hit_g->size() ;
    auto vtx = ev_vtx->at(0); 

    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };
    auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,2);//plane);
         
    auto vtx_w = vtxWT.w ; // Comes in cm; 
    auto vtx_t = vtxWT.t + 800. * geomH->TimeToCm() ; // 800 for single particle files
    //std::cout<<"Vertex info :" <<vtx_w<<", "<<vtx_t<<std::endl ;

    float rad = 0. ;
    float offset = 0.;

    std::vector<float> x(50,0);
    std::vector<float> y(50,0);

    for(int j = 0; j < rad_its; j++){
       //int j = 11;

        rad = float ( (j+1) * 5. ) ;
        _hits_in_rad = 0.;
        _hits_in_rad_g = 0.;

        _hits_p0_in_rad = 0.;
        _hits_p0_in_rad_g = 0.;
        _hits_p1_in_rad = 0.;
        _hits_p1_in_rad_g = 0.;
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
	
	//int h0 = 0;
        //// Count hits in radius
        //for(auto const & h : *ev_hit){

        //    //std::cout<<h.WireID().Wire<<",";

	//    if( h.WireID().Plane != 2 ) continue; 
	//    h0++;

        //    ::cv::Point h_pt(h.WireID().Wire * geomH->WireToCm(),
	//                     h.PeakTime() * geomH->TimeToCm()) ;
        //    auto inside = ::cv::pointPolygonTest(contour,h_pt,false); 
        //   
        //    if(inside  >= 0) _hits_in_rad ++ ;
        //    }

        for(auto const & h : *ev_hit_g){

            ::cv::Point h_pt(h.WireID().Wire * geomH->WireToCm(),
	                     h.PeakTime() * geomH->TimeToCm()) ;
            auto inside = ::cv::pointPolygonTest(contour,h_pt,false); 
        
	    if( h.WireID().Plane == 0 ) {

              if(inside  >= 0) _hits_p0_in_rad_g ++ ;
	      if(inside >= 0 && h.GoodnessOfFit() >= 0 ) _hits_p0_in_rad++;
	      continue;
	    }

	    if( h.WireID().Plane == 1 ) {
              if(inside  >= 0) _hits_p1_in_rad_g ++ ;
	      if(inside >= 0 && h.GoodnessOfFit() >= 0 ) _hits_p1_in_rad++;
	      continue;
	    }

	    if( h.WireID().Plane == 2 ) {

              if(inside  >= 0) _hits_in_rad_g ++ ;
	      if(inside >= 0 && h.GoodnessOfFit() >= 0 ) _hits_in_rad++;
            }
	  }

       // if ( _hits_in_rad_g < 30 ){
       //   bool p1_bigger =  _hits_p1_in_rad_g > _hits_p0_in_rad_g ;

       //   if( p1_bigger ){
       //     _hits_in_rad_g = _hits_p1_in_rad_g ;
       //     _hits_in_rad = _hits_p1_in_rad ;
       //   }
       //   else{
       //     _hits_in_rad_g = _hits_p0_in_rad_g ;
       //     _hits_in_rad = _hits_p0_in_rad ;
       //   }
       // }

        _radii.emplace_back(rad);
	_density.emplace_back(_hits_in_rad / (M_PI * rad * rad )) ;
	_gaus_hits_per_r.emplace_back(_hits_in_rad_g);
	_shr_hits_per_r.emplace_back(_hits_in_rad);

        if( _hits_in_rad_g == 0 )
	  _hits_per_r.emplace_back(0.);
        else {
	  _hits_per_r.emplace_back(float(_hits_in_rad)/_hits_in_rad_g);
	  //std::cout<<"At rad "<<rad<<", ratio: "<<float(_hits_in_rad)/_hits_in_rad_g<<std::endl; // gaus and hit02 ratios: "<<_hits_in_rad_g <<", "<<_hits_in_rad <<std::endl ;
	   //if( float(_hits_in_rad) / _hits_in_rad_g > 0.24)
	   //  _keep++;

	 //if ( _hits_per_r.at(_hits_per_r.size()-1) <= 0.025 && rad == 60 ){ //&& _hits_in_rad_g > 15){

	 if ( _hits_per_r.at(_hits_per_r.size()-1) > 0.975 && rad == 60 && _hits_in_rad_g > 15){
	  
           std::cout<<"\nHIGH RATIO EVENT! "<<_event<<", sub, run, event: "<<storage->subrun_id()<<", "<<storage->run_id()<<", "<<storage->event_id()<<std::endl ;
	   std::cout<<"RATIO : " <<_hits_per_r.at(_hits_per_r.size()-1)<<", and shower hits : "<<_hits_in_rad<<", all: "<<_hits_in_rad_g<<std::endl ;
	 
	   }
	}

      }

    //std::cout<<"Event "<<_event-1<<" has hit_per_r size:  "<<_hits_per_r.size()<<std::endl ;
    //if ( _hits_per_r.at(_hits_per_r.size()-1) > 1. )
     //std::cout<<"RATIO : "<<_hits_per_r.at(_hits_per_r.size()-1)<<std::endl;
    //  std::cout<<"RATIO : "<<_hits_in_rad<<", "<<_hits_in_rad_g<<std::endl;
    //std::cout<<"Event "<<_event-1<<" has hit Dens "<<_hits_per_r.at(8)<<" at radius 45cm "<<std::endl ;
    //std::cout<<"Gauss and shr hits: "<<ev_hit_g->size()<<", "<<ev_hit->size()<<"\n" ;

    _tree->Fill();
  
    return true;
  }

  bool VtxDensity::finalize() {

    //std::cout<<"At 60cm, for 0.24 cut we keep: "<<_keep<<std::endl ;

    if(_fout){
      _fout->cd(); 
      _tree->Write(); 
      }
  
    return true;
  }

}
#endif
