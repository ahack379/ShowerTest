#ifndef LARLITE_SEPTRKSHRNEARVTX_CXX
#define LARLITE_SEPTRKSHRNEARVTX_CXX

#include "SepTrkShrNearVtx.h"
#include "LArUtil/GeometryHelper.h"
#include "Clusterer/Linearity.h"
#include "DataFormat/shower.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/event_ass.h"

namespace larlite {

  bool SepTrkShrNearVtx::initialize() {

    if(!_lin_tree){
      _lin_tree = new TTree("lin_tree","lin_tree");
      _lin_tree->Branch("lin",&_lin,"lin/F");
      _lin_tree->Branch("tll",&_tll,"tll/F");
      _lin_tree->Branch("nhits",&_nhits,"nhits/I");
    }

    _event = 0;
    return true;
  }

  void SepTrkShrNearVtx::Clear(){
    _lin = -999;
    _tll = -999;
    _nhits = -999;
    }
  
  bool SepTrkShrNearVtx::analyze(storage_manager* storage) {

    std::cout<<"\n\nEvent : "<<_event <<std::endl ;
    _event ++ ;
 
    auto geomH  = larutil::GeometryHelper::GetME();

    auto ev_mctrk = storage->get_data<event_mctrack>("mcreco"); 
    auto ev_mcshr = storage->get_data<event_mcshower>("mcreco"); 

   // Want to ignore longest track associated to vertex-- Get handle to association 
    auto ev_ass = storage->get_data<larlite::event_ass>("showerreco");
    auto ev_ass_c = storage->get_data<larlite::event_ass>("ImageClusterHit");

    auto ev_clus = storage->get_data<event_cluster>("ImageClusterHit");
    if ( !ev_clus || !ev_clus->size() ) return false;

    auto ev_hit = storage->get_data<event_hit>("hit02");
    if ( !ev_hit || !ev_hit->size() ) return false;

    auto ev_shr = storage->get_data<event_shower>("showerreco");
    if ( !ev_shr || !ev_shr->size() ) return false;

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");
    if ( !ev_vtx || !ev_vtx->size() ) return false;

    // Get association to shr => clus 
    auto const& ass_clus_v = ev_ass->association(ev_shr->id(), ev_clus->id());
    if (!ass_clus_v.size()) {
      std::cout << "No ass from shower -> clus! " << std::endl;
      return false;
      }   

    // Get association to shr => clus 
    auto const& ass_hit_v = ev_ass_c->association(ev_clus->id(), ev_hit->id());
    if (!ass_hit_v.size()) {
      std::cout << "No ass from clus -> hit! " << std::endl;
      return false;
      }

    auto vtx = ev_vtx->at(0);
    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };

    std::vector<int> clus_ids ;
    std::vector<int> plane_ids ;
    std::vector<int> shr_ids ;

    for( size_t s = 0; s < ev_shr->size(); s++) { // auto const & s : ev_shr ){
      auto const shr = ev_shr->at(s) ;

      for( size_t i = 0; i < ass_clus_v.at(s).size(); i++){
        auto clus = ev_clus->at(ass_clus_v.at(s).at(i)) ;

        //std::cout<<"Plane: "<<clus.Plane().Plane<<std::endl ;

        auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,clus.Plane().Plane);
        auto vtx_w = vtxWT.w;
        auto vtx_t = vtxWT.t + 800 * geomH->TimeToCm() ;

        auto shr_st_2d = geomH->Point_3Dto2D(shr.ShowerStart(),clus.Plane().Plane); 

        auto dist_to_vtx = sqrt( pow( clus.StartWire() * geomH->WireToCm() - vtx_w,2) + 
                                 pow( clus.StartTick() * geomH->TimeToCm() - vtx_t,2) );
        
        //std::cout<<"Time : "<<clus.StartTick()<<", "<<clus.StartTick()*geomH->TimeToCm()<<", "<<vtx_t<<std::endl ;
        // If the distance to the vertex is too big, don't consider this "shower" in our sample
        if ( dist_to_vtx > 3.5 ) continue; //std::cout<<"Got one! "<<std::endl ;

        clus_ids.emplace_back(ass_clus_v.at(s).at(i));
        plane_ids.emplace_back(clus.Plane().Plane);
        shr_ids.emplace_back(s);
       } 
     }

     std::vector<int> mcshr_ids ;
     std::vector<int> mctrk_ids;
     for ( size_t p = 0; p < plane_ids.size(); p++){ 

       auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,plane_ids[p]); //it->first);
       auto vtx_w = vtxWT.w;
       auto vtx_t = vtxWT.t + 800 * geomH->TimeToCm() ;

       std::cout<<"In Plane : "<<plane_ids[p]<<std::endl ;
       std::cout<<"VERTEx : "<<vtx_w/geomH->WireToCm()<<", "<<vtx_t/geomH->TimeToCm()<<std::endl;
       for ( size_t k = 0; k < ev_mcshr->size(); k++){ 

         auto traj = ev_mcshr->at(k).DetProfile();
         auto p_st_2d = geomH->Point_3Dto2D(traj.X(),traj.Y(),traj.Z(),plane_ids[p]); //it->first); 

	 auto dist_to_vtx = sqrt( pow( p_st_2d.w - vtx_w,2) + 
                                  pow( p_st_2d.t - vtx_t,2) );
	

	// if (traj.E() > 30 && fabs(vtx_w - p_st_2d.w)*geomH->WireToCm() < 5) { // && p == 0) 
	//   std::cout<<"MCS St: "<<p_st_2d.w/geomH->WireToCm()<<", "<<p_st_2d.t/geomH->TimeToCm() + 800<<std::endl;
	// }
        
         if( dist_to_vtx > 3.5 ) continue; 
	   std::cout<<"Dist : "<<dist_to_vtx <<std::endl ;

	 //std::cout<<"FOUDN ONE "<<std::endl ;

         mcshr_ids.emplace_back(k) ;
       }

     for ( size_t k = 0; k < ev_mctrk->size(); k++){ 

         auto traj = ev_mctrk->at(k).Start();
         auto p_st_2d = geomH->Point_3Dto2D(traj.X(),traj.Y(),traj.Z(),plane_ids[p]); //it->first); 

	 auto dist_to_vtx = sqrt( pow( p_st_2d.w - vtx_w,2) + 
                                  pow( p_st_2d.t - vtx_t,2) );
         if( dist_to_vtx > 3.5 ) continue; 

         mctrk_ids.emplace_back(k) ;
       }
     }

     std::cout<<"Size of reco and mc : "<<shr_ids.size()<<", "<<mcshr_ids.size()<<std::endl; 

     // Create map to hold largest scores; maximizing dot product to match
     // mc to reco shower
     std::multimap<float,std::pair<int,int>> score_to_match ;

     int prev_id = -1 ;

     for( int ii = 0 ; ii < shr_ids.size(); ii++ ){

       //Store shower IDs with cluster pair projections; want to only consider 1 shower
       if( shr_ids[ii] == prev_id ) continue; 
     
       auto shr_ii = ev_shr->at(shr_ids[ii]);
       auto reco_st = shr_ii.ShowerStart();
       auto mag_ii = sqrt(pow(reco_st.X(),2) + pow(reco_st.Y(),2) + pow(reco_st.Z(),2) );
       prev_id = shr_ids[ii];

       for( int jj = 0 ; jj < mcshr_ids.size()+mctrk_ids.size(); jj++ ){
         // Store shower ids for each cluster in each plane; so will have shower duplicates
         if( jj < mcshr_ids.size() ){

           auto mcshr_jj = ev_mcshr->at(mcshr_ids[jj]); 
	   auto mc_st = mcshr_jj.StartDir() ;

	   auto mag_jj = sqrt(pow(mc_st.X(),2) + pow(mc_st.Y(),2) + pow(mc_st.Z(),2) );
	   auto dot = (mc_st.X()*reco_st.X() + mc_st.Y()*reco_st.Y() + mc_st.Z() * reco_st.Z())/mag_jj/mag_ii ;

//	   std::cout<<"Dir mc shower: "<<mcshr_jj.StartDir().X()<<", "<<mcshr_jj.StartDir().Y()<<", "<<mcshr_jj.StartDir().Z()<<std::endl ;
//	   std::cout<<"Dir mc shower: "<<shr_ii.Direction().X()<<", "<<shr_ii.Direction().Y()<<", "<<shr_ii.Direction().Z()<<std::endl ;
	   //std::cout<<"SHOWER DOt is: "<<dot<<std::endl ;
	   std::pair<int,int> recomc = { ii, jj }; 
	   score_to_match.emplace(1./dot,recomc) ;//(std::make_pair<ii,jj>));
	 }
         else{
	   auto jj_adj = jj - mcshr_ids.size() ;
           auto mctrk_jj = ev_mctrk->at(mctrk_ids[jj_adj]); 
	   auto mc_st = mctrk_jj.Start() ;

	   auto mag_jj = sqrt(pow(mc_st.X(),2) + pow(mc_st.Y(),2) + pow(mc_st.Z(),2) );
	   auto dot = (mc_st.X()*reco_st.X() + mc_st.Y()*reco_st.Y() + mc_st.Z() * reco_st.Z())/mag_jj/mag_ii ;

	   std::cout<<"TRACK DOt is: "<<dot<<std::endl ;
	   std::pair<int,int> recomc = { ii, jj }; 
	   score_to_match.emplace(1./dot,recomc) ;//(std::make_pair<ii,jj>));
          } 
         }
       }

       //for( auto const & i : score_to_match)
       //  std::cout<<"Score: "<<i.first<<", reco + mc id: "<<i.second.first<<", "<<i.second.second<<std::endl ;

       // Now we have a map of dot scores with reco + mc trk IDS; 

    for ( size_t j = 0; j < clus_ids.size(); j++ ){

      Clear();

      auto shr_i = ev_shr->at(shr_ids[j]);

      std::vector<double> wire_v ; 
      std::vector<double> time_v ; 

      _nhits = ass_hit_v.at(clus_ids[j]).size() ;

      for ( size_t k = 0; k < _nhits; k++ ){

         auto h = ev_hit->at(k) ;
 
         wire_v.emplace_back(h.WireID().Wire);
         time_v.emplace_back(h.PeakTime());
         }   

      ::Linearity linear(wire_v,time_v);
      linear._r = 10 ;
   
      linear.linearity(wire_v,time_v) ;
      linear.local_linearity(wire_v,time_v) ;
    
      _lin = linear._lin ;
      _tll = linear._local_lin_truncated_avg;

      _lin_tree->Fill();
      }

 

    return true;
  }

  bool SepTrkShrNearVtx::finalize() {

    if(_fout) { _fout->cd(); _lin_tree->Write(); }
  
    return true;
  }

}
#endif
