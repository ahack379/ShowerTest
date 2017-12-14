#ifndef LARLITE_TESTSlimmedMCPART_CXX
#define LARLITE_TESTSlimmedMCPART_CXX

#include "TestSlimmedMCPart.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/hit.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/cluster.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"

namespace larlite {

  bool TestSlimmedMCPart::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //

    return true;
  }
  
  bool TestSlimmedMCPart::analyze(storage_manager* storage) {

    std::cout<<"New event : "<<std::endl ;

    auto ev_ass = storage->get_data<larlite::event_ass>("gaushitTruthMatch");

    std::cout<<"Hit association tree size: "<<ev_ass->size()<<std::endl ;
     if ( !ev_ass || ev_ass->size() == 0 ) {
       std::cout << "No such association! " << std::endl;
       return false;
     }

    auto const& ass_hit_v = ev_ass->association(0);
    auto const& ass_mcpart_v = ev_ass->association(1);

    auto ev_mcshower = storage->get_data<event_mcshower>("mcreco");
    auto ev_mct = storage->get_data<event_mctrack>("mcreco");

    auto ev_hit = storage->get_data<event_hit>("gaushit");

    auto ev_mcpart = storage->get_data<event_mcpart>("largeant");

    // grab associations by product id
    auto const& ass_hit_mcpart_v = ev_ass->association(ev_hit->id(),ev_mcpart->id());

    // Retrieve cluster data product (output)
    auto ev_mccluster = storage->get_data<event_cluster>("mccluster_slim");
    auto cluster_ass_v = storage->get_data<event_ass>(ev_mccluster->name());

    storage->set_id(ev_hit->run(),ev_hit->subrun(),ev_hit->event_id());

    //std::map<int,std::vector<int>> mcpart_to_hit_map ;
    std::vector<std::vector<unsigned int> > cluster_hit_v;
    cluster_hit_v.resize( 3 *( ev_mcshower->size() + ev_mct->size()) );
    std::vector<larlite::geo::View_t> cluster_plane_v(cluster_hit_v.size(), larlite::geo::View_t::kUnknown);

    // mcparticle (TrackId) to mcshower map
    std::map<int,int> mcp_to_mcs ;

    std::map<int,int> mcp_to_mcst ;

    for (int jj=0; jj < ev_mcshower->size(); jj++){

      auto mcs = ev_mcshower->at(jj); 
      if (mcs.DetProfile().E() != 0 ){
        for ( auto const & p : mcs.DaughterTrackID() ){

          //std::cout<<"Energy : "<<mcs.DetProfile().E() <<" trackid: "<<p<<std::endl ;
          mcp_to_mcst[p] = jj ;
        }
      }

    }

    // mcparticle to mctrack map
    std::map<int,int> mcp_to_mct ;

    for (int jj=0; jj < ev_mct->size(); jj++){
      
      auto mct = ev_mct->at(jj); 
      if ( mcp_to_mcst.find(mct.TrackID() ) != mcp_to_mcst.end() ){
         std::cout<<"So...we can have overlap in trackIDs? "<<mct.TrackID()<<std::endl ;
         continue;
      }
      if (mct.size() < 2 or (mct.Origin() != 1 and mct.Origin() != 2) ){ continue; }

      mcp_to_mcst[mct.TrackID()] = (jj + ev_mcshower->size() ) ;
    }


    // loop over hits
    for (size_t kk=0; kk < ass_hit_mcpart_v.size(); kk++) {

      int hit_idx = kk ;
      auto const& hit    = ev_hit->at(hit_idx);
      auto const& pl     = hit.View();
      auto const& mcpart_id = ass_hit_mcpart_v.at(kk);

      if ( mcpart_id.size() == 0 ) continue; 
      int mcp_trackid = mcpart_id.at(0);
     
      if ( mcp_to_mcst.find(mcp_trackid) != mcp_to_mcst.end() ){
        int idx = mcp_to_mcst[mcp_trackid] ;
        cluster_hit_v[ pl * (ev_mcshower->size() + ev_mct->size() ) + idx  ].push_back( hit_idx );
        cluster_plane_v[ pl *(ev_mcshower->size() + ev_mct->size() ) + idx ] = pl;
	//std::cout<<"MC shower trackID: "<<mcp.TrackId()<<", "<<ev_mcshower->at(mcp_to_mcs[mcp.TrackId()]).TrackID()<<std::endl ;
      }
      //else if ( mcp_to_mct.find(mcp_trackid) != mcp_to_mct.end() ){
      //  int idx = mcp_to_mct[mcp_trackid] ;
      //  cluster_hit_v[ pl * ( ev_mcshower->size()) + ev_mcshower->size() + idx  ].push_back( hit_idx );
      //  cluster_plane_v[ pl *( ev_mcshower->size()) + ev_mcshower->size() + idx ] = pl;
      //}
      else continue;
      //std::cout<<"Hit info : "<<ev_hit->at(kk).WireID().Plane <<std::endl ;
      //std::cout<<"Particle : "<<mcp.PdgCode()<<", "<<mcp.TrackId()<<", "<<std::endl;
    }

    std::cout<<" CLusters: "<<cluster_hit_v.size() <<std::endl;

    std::vector<std::vector<unsigned int> > cluster_hit_ass_v;
    for (size_t idx=0; idx < cluster_hit_v.size(); idx++){
      if (cluster_hit_v[idx].size() < 4) 
        continue;

      // create a new cluster
      cluster clus;
      clus.set_n_hits(cluster_hit_v[idx].size());
      clus.set_view(cluster_plane_v[idx]);
      ev_mccluster->push_back(clus);
      cluster_hit_ass_v.push_back(cluster_hit_v[idx]);

    }// for all clusters created

    // now save the associations for the cluster
    cluster_ass_v->set_association(ev_mccluster->id(),product_id(data::kHit,ev_hit->name()),
                                   cluster_hit_ass_v);


    return true;
  }

  bool TestSlimmedMCPart::finalize() {

    return true;
  }

}
#endif
