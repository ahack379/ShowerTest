#ifndef LARLITE_FILTERHITRATIO_CXX
#define LARLITE_FILTERHITRATIO_CXX

#include "FilterHitRatio.h"
#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"
#include "LArUtil/GeometryHelper.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <math.h>

namespace larlite {

  bool FilterHitRatio::initialize() {

   _event = 0 ;
   _hits_in_rad_g = 0;
   _hits_in_rad = 0;

   if( !_tree ){
     _tree = new TTree("tree","tree");
     _tree->Branch("hit_ratio","std::vector<float>",&_hit_ratio); 
    }   

    return true;
  }
  
  bool FilterHitRatio::analyze(storage_manager* storage) {


    auto const& geomH = ::larutil::GeometryHelper::GetME();

    std::cout<<"\nNew event! "<<_event<<std::endl ;
    _event++;
    
    auto ev_hit = storage->get_data<event_hit>("shrhits");
    if ( !ev_hit || !ev_hit->size() ) return false;

    auto ev_hit_g = storage->get_data<event_hit>("gaushit");
    if ( !ev_hit_g || !ev_hit_g->size() ) return false; 

    _hit_ratio.clear();

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");

    if ( !ev_vtx || !ev_vtx->size() ) return false; 

    auto vtx = ev_vtx->at(0); 

    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };
    auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,2);//plane);
    
    auto vtx_w = vtxWT.w; // / geomH->WireToCm();
    auto vtx_t = vtxWT.t + 800 * geomH->TimeToCm() ; // geomH->TimeToCm() + 800; //_time_offset ;

    std::vector<::cv::Point> contour;

    std::vector<float> x(50,0);
    std::vector<float> y(50,0);

    // Build circle around radius
    for(int i = 0; i < x.size(); i++){

       x[i] = vtx_w + _radius * cos(M_PI * 2. * i / x.size());
       y[i] = vtx_t + _radius * sin(M_PI * 2. * i / x.size());

       ::cv::Point pt(x[i],y[i]) ;
       contour.emplace_back(pt);
       }

     std::cout<<"Gauss and shr hits: "<<ev_hit_g->size()<<", "<<ev_hit->size()<<std::endl ;

     // Count shower hits in radius
     for(int i = 0; i < ev_hit_g->size(); i++){

         auto const & h_g = ev_hit_g->at(i); 
      
         ::cv::Point h_pt_g(h_g.WireID().Wire * geomH->WireToCm(),
                          h_g.PeakTime() * geomH->TimeToCm()) ;
         auto inside_g = ::cv::pointPolygonTest(contour,h_pt_g,false);

         if(inside_g >= 0) _hits_in_rad_g++ ;
         
         if( i < ev_hit->size() ){

             auto const & h = ev_hit->at(i);
             
             ::cv::Point h_pt(h.WireID().Wire * geomH->WireToCm(),
                              h.PeakTime() * geomH->TimeToCm()) ;
             auto inside = ::cv::pointPolygonTest(contour,h_pt,false);

             if(inside  >= 0) _hits_in_rad ++ ;
             }
       }

    auto ev_truth = storage->get_data<event_mctruth>("generator");
    auto & truth = ev_truth->at(0);
    auto & nu  = truth.GetNeutrino();
    
    auto const & intType = nu.InteractionType();
    if( IsSignal(intType) ) _total_sig_events++ ;

     if( _hits_in_rad_g == 0 ){
       //_hit_ratio.emplace_back(0.);
       return false ;
        }
     else{
       float ratio = float( _hits_in_rad ) / _hits_in_rad_g ;
       //_hit_ratio.emplace_back( ratio );
        
       if( ratio < _ratio_cut ){
          return false; 
         }
       else 
          if( IsSignal(intType) )
            _good_tags ++ ;
          else
            _bad_tags ++ ;
          return true; 
        }

     //std::cout<<"Hits in rad "<<rad<<": "<<float(_hits_in_rad)/_hits_in_rad_g<<std::endl; 

    return true;
  }

  bool FilterHitRatio::finalize() {

    //if(_fout) { _fout->cd(); _tree->Write(); }

    std::cout<<float(_good_tags)/_total_sig_events*100<<"\% of signal events in Selection 2 pass this filter"<<std::endl ;
    std::cout<<float(_bad_tags)/(_good_tags + _bad_tags)*100<<"\% of our final sample is shit we don't want"<<std::endl ;
  
    return true;
  }


// t is interaction type
bool FilterHitRatio::IsSignal(int t){

    if( (t == 1004 || t == 1011 || t == 1080 || t == 1086 || t == 1090) )  
      return true ;
    else 
      return false ;
  }

}
#endif
