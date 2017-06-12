#ifndef LARLITE_SHOWERCOMPLETENESS_CXX
#define LARLITE_SHOWERCOMPLETENESS_CXX

#include "ShowerCompleteness.h"
#include "DataFormat/cluster.h"
#include "DataFormat/shower.h"
#include "DataFormat/hit.h"

namespace larlite {

  bool ShowerCompleteness::initialize() {

    // Fill completeness / purity tree
    if ( !_tree ){
      _tree = new TTree("tree","tree");
      _tree->Branch("_purity_v","std::vector<std::vector<float>>",&_purity_v);
      _tree->Branch("_complete_v","std::vector<std::vector<float>>",&_complete_v);
    }

   _p0 = 0;
   _p1 = 0; 
   _p2 = 0;
   _tot_clus = 0;
   _event = -1;

   _purity_v.resize(3) ;
   _complete_v.resize(3) ;

    return true;
  }
  
  bool ShowerCompleteness::analyze(storage_manager* storage) {

    _event++;
    //std::cout<<"EVENT "<<std::endl ;

    _mc_hit_map.clear();

    auto ev_ass = storage->get_data<larlite::event_ass>("mccluster");
    auto const& ass_keys = ev_ass->association_keys();

    if ( ass_keys.size() == 0 ) return false; 

    larlite::event_cluster *ev_mcclus = nullptr;
    auto ass_hit_clus_v = storage->find_one_ass( ass_keys[0].first, ev_mcclus, ev_ass->name() );

    larlite::event_hit *ev_mchit = nullptr;
    auto ass_mcclus_v = storage->find_one_ass( ass_keys[0].second, ev_mchit, ev_ass->name() );

    auto ev_hit = storage->get_data<larlite::event_hit>("gaushit");

    if (ass_hit_clus_v.size() == 0){
      std::cout << "No ass! exit" << std::endl;
      return false;
    }   
    if (ass_mcclus_v.size() == 0){
      std::cout << "No ass! exit" << std::endl;
      return false;
    }   

    if (!ev_hit || ev_hit->size() == 0){
      std::cout << "No ass! exit" << std::endl;
      return false;
    }   

   for( int i = 0; i < ass_mcclus_v.size(); i ++ ){
     for ( int j = 0; j < ass_mcclus_v[i].size(); j++ ){

       auto hid = ass_mcclus_v[i][j];
       _mc_hit_map[hid] = i ; 
     }

   }

   // Now we have the mc hit info in a map compute purity and completeness

   // Get the association from cluster -> hit
   auto const & ev_clus = storage->get_data<event_cluster>("ImageClusterHit");
   auto const & ev_ass_c = storage->get_data<larlite::event_ass>("ImageClusterHit");
   auto const & ass_imageclus_v = ev_ass_c->association(ev_clus->id(), ev_hit->id());

   // Get the association from shower -> cluster
   auto const & ev_s = storage->get_data<event_shower>("pi0_candidate_showers"); 
   auto ev_ass_s = storage->get_data<larlite::event_ass>("pi0_candidate_showers");
   auto const& ass_showerreco_v = ev_ass_s->association(ev_s->id(), ev_clus->id());

   if ( !ev_clus || ev_clus->size() == 0 || !ev_ass_c || ev_ass_c->size() == 0 ){
     std::cout<<"No cluster..." <<std::endl; return false;
     }   
   //if ( !ass_imageclus_v || ass_imageclus_v.size() == 0 ){ std::cout<<"No imageclus ass..." <<std::endl; return false; }   

   if( !ev_s || !ev_s->size() ) { 
     std::cout<<"No shower..." <<std::endl;
     return false;
     }   

   // for each reco cluster, find the origin of all hits and calc purity/completeness 
   std::vector<int> pur_ctr_v ;

   std::cout<<std::endl ;

   // Loop over showers. I think.
   for (size_t i = 0; i < ass_showerreco_v.size(); i++ ){
     
     _tot_clus += ass_showerreco_v.at(i).size() ;
     //std::cout<<"Number of clusters assoc to shower: "<<ass_showerreco_v.at(i).size()<<std::endl ;
     // Loop over clusters associated to this shower. Again, I think.
     for (size_t j = 0; j < ass_showerreco_v.at(i).size(); j++ ){

         auto clus_id = ass_showerreco_v.at(i).at(j); 
         auto iclus = ev_clus->at(clus_id);
     
         int plane = iclus.Plane().Plane ;

         pur_ctr_v.resize(ass_mcclus_v.size(),0) ;

         int max_hits = -1;
         int max_cid = -1 ;

         //auto const & ass_clus_hit = ass_imageclus_v.at(clus_id) ;

         for ( int k = 0; k < ass_imageclus_v.at(clus_id).size(); k++ ){

           //auto hid = ass_imageclus_v[j][k] ;
           auto hid = ass_imageclus_v.at(clus_id).at(k) ; 
           auto h = ev_hit->at(hid);
           
           if ( _mc_hit_map.find(hid) != _mc_hit_map.end() ){
             pur_ctr_v[_mc_hit_map[hid]]++ ; 
             if( pur_ctr_v[ _mc_hit_map[hid]] > max_hits ){
               max_hits = pur_ctr_v[_mc_hit_map[hid]];
               max_cid = _mc_hit_map[hid] ;
             }
           }
         }

	 float purity = 0 ;
	 float complete = 0;

         if ( max_cid != -1 ){

           auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
           auto tot_reco_hits = ass_imageclus_v[clus_id].size();
           
           purity   = float(max_hits) / tot_reco_hits ;
           complete = float(max_hits) / tot_mc_hits ;


	   if ( purity > 1 ){ 
	     std::cout<<"****************************************************************"<<std::endl;
	     }
	     //std::cout<<"CID: "<<max_cid <<std::endl ;
	     std::cout<<"Event: "<<_event<<", plane: "<<plane<<std::endl ;
	     std::cout<<"Purtiy : "<<purity<<", "<<complete<<std::endl ;
	     std::cout<<"Max hits, Reco Tot, MC Tot : "<<max_hits<<", "<<tot_reco_hits<<", "<<tot_mc_hits<<std::endl ;
	  //}
	}

         _purity_v[plane].emplace_back(purity) ;
         _complete_v[plane].emplace_back(complete) ;

	 if( plane == 0 )
	   _p0++ ;
	 else if ( plane == 1) 
	   _p1++ ;
	 else if ( plane == 2)
	   _p2++ ;

	//_tot_clus ++ ;
      }
   }

   //for ( int j = 0; j < ass_imageclus_v.size(); j++ ){

   //  pur_ctr_v.resize(ass_mcclus_v.size(),0) ;

   //  int max_hits = -1;
   //  int max_cid = -1 ;
   //  int plane = -1 ;

   //  for ( int k = 0; k < ass_imageclus_v[j].size(); k++ ){

   //    auto hid = ass_imageclus_v[j][k] ;
   //    auto h = ev_hit->at(hid); 
   //    plane = h.WireID().Plane ;
   //    
   //    if ( _mc_hit_map.find(hid) != _mc_hit_map.end() ){
   //      pur_ctr_v[_mc_hit_map[hid]]++ ; 
   //      if( pur_ctr_v[ _mc_hit_map[hid]] > max_hits ){
   //        max_hits = pur_ctr_v[_mc_hit_map[hid]];
   //        max_cid = _mc_hit_map[hid] ;
   //      }
   //    }
   //  }
   //  auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
   //  auto tot_reco_hits = ass_imageclus_v[j].size();
   //  
   //  auto purity   = float(max_hits) / tot_reco_hits ;
   //  auto complete = float(max_hits) / tot_mc_hits ;

   //  
   //  _purity_v[plane].emplace_back(purity) ;
   //  _complete_v[plane].emplace_back(complete) ;
   //}

   //
   // 0) Loop through MC pi0 showers that live in the file 
   // 1) For each MC shower, fill an mc hit map per plane 
   // 2) Next find the corresponding reco shower to the MC by taking dot product
   // 3) Compare the hits in the reco shower to the hits filled in the map 
   // 4) For each plane and for each cluster in each plane, compute the purity + completeness 
   //    of the cluster

    return true;
  }

  bool ShowerCompleteness::finalize() {

    _tree->Fill();

    std::cout<<" Total fills: "<<_tot_clus<<std::endl ;
    std::cout<<"By plane : "<<_p0<<", "<<_p1<<", "<<_p2<<std::endl ;

    if( _fout ){
      _fout->cd();
      _tree->Write();
    }
  
    return true;
  }

}
#endif
