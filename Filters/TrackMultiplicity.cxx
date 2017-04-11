#ifndef LARLITE_TRACKMULTIPLICITY_CXX
#define LARLITE_TRACKMULTIPLICITY_CXX

#include "TrackMultiplicity.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"

namespace larlite {

  bool TrackMultiplicity::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool TrackMultiplicity::analyze(storage_manager* storage) {
  
     auto ev_trk = storage->get_data<event_track>("pandoraNu");
    if ( !ev_trk || !ev_trk->size() ) {std::cout<<"No Track!" <<std::endl ; return false; }

    auto ev_recov= storage->get_data<event_vertex>("numuCC_vertex");
    if(!ev_recov || !ev_recov->size() ) {
      std::cout<<"Event has no recovertex info "<<std::endl;
      return false;
      }

    auto vtx = ev_recov->at(0); 
    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };

    //Map of lengths -> track id
    std::multimap<float,int> trk_map ;
    float min_dist = 10000;

    // Find closest + longest pandoraNu track to vertex
    for ( size_t ti = 0; ti < ev_trk->size(); ti++ ) { 

      auto t_vtx = ev_trk->at(ti).Vertex() ;
      auto t_end = ev_trk->at(ti).End() ;
    
      float dist_st = sqrt( pow(t_vtx.X() - vtxXYZ[0],2) + 
                            pow(t_vtx.Y() - vtxXYZ[1],2) + 
                            pow(t_vtx.Z() - vtxXYZ[2],2) );  

      float dist_end = sqrt( pow(t_end.X() - vtxXYZ[0],2) + 
                             pow(t_end.Y() - vtxXYZ[1],2) + 
                             pow(t_end.Z() - vtxXYZ[2],2) );  
       if ( dist_st < 3 || dist_end < 3 ){
          float len = ev_trk->at(ti).Length();
          trk_map.emplace(1./len,ti);
          min_dist = dist_st < dist_end ? dist_st : dist_end ; 
        }   
     }   

    if( !trk_map.size() ) return false;

    if( ev_trk->at(trk_map.begin()->second).Length() < 15 )
        return false;
 
    return true;
  }

  bool TrackMultiplicity::finalize() {

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
    // if(_fout) { _fout->cd(); h1->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    return true;
  }

}
#endif
