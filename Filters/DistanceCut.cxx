#ifndef LARLITE_DISTANCECUT_CXX
#define LARLITE_DISTANCECUT_CXX

#include "DistanceCut.h"
#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "LArUtil/GeometryHelper.h"
#include <math.h>

namespace larlite {

  bool DistanceCut::initialize() {

   if( !_dist_cut_tree ){
     _dist_cut_tree = new TTree("dist_cut_tree","dist_cut_tree"); 
     _dist_cut_tree->Branch("hits_in_rad",&_hits_in_rad,"hits_in_rad/F"); 
     _dist_cut_tree->Branch("hits_in_rad_g",&_hits_in_rad_g,"hits_in_rad_g/F"); 
     _dist_cut_tree->Branch("hits_tot",&_hits_tot,"hits_tot/F"); 
     _dist_cut_tree->Branch("event",&_event,"event/I"); 
     _dist_cut_tree->Branch("hits_per_r","std::vector<float>",&_hits_per_r); 
    }

    _event = 0;

    return true;
  }
  
  bool DistanceCut::analyze(storage_manager* storage) {

    int rad_its = 1;
    _hits_tot = 0.;
    _hits_per_r.clear();
    _hits_per_r.reserve(rad_its);

    auto const& geomH = ::larutil::GeometryHelper::GetME();

    //std::cout<<"\nNew event! "<<_event<<std::endl ;
    _event++;

    auto ev_hit_g = storage->get_data<event_hit>("gaushit");
    if ( !ev_hit_g || !ev_hit_g->size() ) {std::cout<<"Returning, no hits..."<<std::endl ; return false; }

    auto ev_vtx = storage->get_data<event_vertex>("mcvertex"); //numuCC_vertex");
    if ( !ev_vtx || !ev_vtx->size() ) {std::cout<<"Returning, no vertex..."<<std::endl ; return false; }

    auto vtx = ev_vtx->at(0); 

    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };
    auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,2);//plane);
         
    auto vtx_w = vtxWT.w ; // Comes in cm; 
    auto vtx_t = vtxWT.t + 800. * geomH->TimeToCm() ; // 800 for single particle files
 
    float ratio = 0.;
    _hits_in_rad = 0;
    _hits_in_rad_g = 0;

    for(auto const & h : *ev_hit_g){

      if( h.WireID().Plane != 2 ) continue; 

        auto w = h.WireID().Wire * geomH->WireToCm();
        auto t = h.PeakTime() * geomH->TimeToCm() ;

	auto dist = sqrt( pow(w - vtx_w,2) + pow(t - vtx_t,2) );

        if(dist <= _radius ) _hits_in_rad_g ++ ;

        if(dist <= _radius && h.GoodnessOfFit() >= 0) _hits_in_rad++ ;

	//std::cout<<"Dist and rad: "<<dist<<", "<<_radius <<std::endl ;
     }
        
     if( _hits_in_rad_g < 15 ) 
	ratio = 0;
     else 
	ratio = float(_hits_in_rad)/_hits_in_rad_g;

     //std::cout<<"Event "<<_event-1<<" has Gaus and shr hits: "<<_hits_in_rad_g <<", "<<_hits_in_rad <<std::endl ;

    if( ratio < _ratio_cut ) return false ;

    _hits_tot = ev_hit_g->size() ;
    _dist_cut_tree->Fill();
  
    return true;
  }

  bool DistanceCut::finalize() {

    std::cout<<"\n\n****** "<<_dist_cut_tree->GetEntries()<<" entries have passed the distance cut ******\n"<<std::endl ;

    if(_fout){
      _fout->cd(); 
      _dist_cut_tree->Write(); 
      }
  
    return true;
  }

}
#endif
