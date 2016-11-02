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
     _tree->Branch("hits_in_rad",&_hits_in_rad,"hits_in_rad/F"); 
     _tree->Branch("hits_in_rad_g",&_hits_in_rad_g,"hits_in_rad_g/F"); 
     _tree->Branch("hits_tot",&_hits_tot,"hits_tot/F"); 
     _tree->Branch("radii","std::vector<float>",&_radii); 
     _tree->Branch("density","std::vector<float>",&_density); 
     _tree->Branch("event",&_event,"event/I"); 
     _tree->Branch("hits_per_r","std::vector<float>",&_hits_per_r); 
    }

    _event = 0;

    return true;
  }
  
  bool VtxDensity::analyze(storage_manager* storage) {

    int rad_its = 15;
    _hits_tot = 0.;

    _radii.clear();
    _density.clear();
    _hits_per_r.clear();

    _radii.reserve(rad_its);
    _density.reserve(rad_its);
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

    float rad = 0. ;
    float offset = 0.;

    std::vector<float> x(50,0);
    std::vector<float> y(50,0);

    std::cout<<"\nvtx coords : " <<vtx_w<<", "<<vtx_t<<", and wire time: "<<vtx_w/geomH->WireToCm()<<", "<<vtx_t/geomH->TimeToCm()<<std::endl ;
    
    for(int j = 0; j < rad_its; j++){

        rad = (j+1) * 5 + offset ;
        _hits_in_rad = 0.;
        _hits_in_rad_g = 0.;

        // Build circle around radius
        for(int i = 0; i < x.size(); i++){

           x[i] = vtx_w + rad * cos(M_PI * 2. * i / x.size());
           y[i] = vtx_t + rad * sin(M_PI * 2. * i / x.size());

	   //std::cout<<"["<<x[i]<<", "<<y[i]<<"]," ;

           ::cv::Point pt(x[i],y[i]) ;
           contour.emplace_back(pt);
           } 

        ::cv::Point pt(contour.at(0).x, contour.at(0).y) ;
        contour.emplace_back(pt);
	

        // Count hits in radius
        for(auto const & h : *ev_hit){

            //std::cout<<h.WireID().Wire<<",";

	    if( h.WireID().Plane != 2 ) continue; 
            //std::cout<<"Hit loc + vtx : ("<<h.WireID().Wire<<","<<h.PeakTime()<<"),("<<vtx_w<<","<<vtx_t<<")"<<std::endl;

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

	std::cout<<"At rad "<<rad<<" gaus and hit02 ratios: "<<_hits_in_rad_g <<", "<<_hits_in_rad <<std::endl ;

        _radii.emplace_back(rad);
	_density.emplace_back(_hits_in_rad / (M_PI * rad * rad )) ;

      if( _use_mcbnb_info ){
        auto ev_truth = storage->get_data<event_mctruth>("generator");
        auto & truth = ev_truth->at(0);
        auto & nu  = truth.GetNeutrino();
        
        auto const & t = nu.InteractionType();
        if( (t == 1004 || t == 1011 || t == 1080 || t == 1086 || t == 1090) ){

          if( _hits_in_rad_g == 0 )
            _hits_per_r.emplace_back(0.);
          else 
            _hits_per_r.emplace_back(float(_hits_in_rad)/_hits_in_rad_g);
            } 
          }
      else{

          if( _hits_in_rad_g == 0 )
	    _hits_per_r.emplace_back(0.);
          else 
	    _hits_per_r.emplace_back(float(_hits_in_rad)/_hits_in_rad_g);
        //std::cout<<"Hits in rad "<<rad<<": "<<float(_hits_in_rad)/2/_hits_in_rad_g<<std::endl; 
        }

      }

    std::cout<<"Event "<<_event-1<<" has hit Dens "<<_hits_per_r.at(8)<<" at radius 45cm "<<std::endl ;
    std::cout<<"Gauss and shr hits: "<<ev_hit_g->size()<<", "<<ev_hit->size()<<"\n" ;



    
    _hits_tot = ev_hit->size() ;
    _tree->Fill();
  
    return true;
  }

  bool VtxDensity::finalize() {

    if(_fout){
      _fout->cd(); 
      _tree->Write(); 
      }
  
    return true;
  }

}
#endif
