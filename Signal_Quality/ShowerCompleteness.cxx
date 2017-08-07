#ifndef LARLITE_SHOWERCOMPLETENESS_CXX
#define LARLITE_SHOWERCOMPLETENESS_CXX

#include "ShowerCompleteness.h"
#include "DataFormat/cluster.h"
#include "DataFormat/shower.h"
#include "DataFormat/hit.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctruth.h"


namespace larlite {

  bool ShowerCompleteness::initialize() {

    // Fill completeness / purity tree
    if ( !_tree ){
      _tree = new TTree("tree","tree");
      _tree->Branch("_purity_v","std::vector<std::vector<float>>",&_purity_v);
      _tree->Branch("_complete_v","std::vector<std::vector<float>>",&_complete_v);
      _tree->Branch("_cw_purity_v","std::vector<std::vector<float>>",&_cw_purity_v);
      _tree->Branch("_cw_complete_v","std::vector<std::vector<float>>",&_cw_complete_v);
    }

    if ( !_pi0 ){
      _pi0 = new TTree("pi0","pi0");
      _pi0->Branch("mass",&_mass,"mass/F");
      _pi0->Branch("min_pur",&_min_pur,"min_pur/F");
      _pi0->Branch("min_com",&_min_com,"min_com/F");
      _pi0->Branch("shr1_dot",&_shr1_dot,"shr1_dot/F");
      _pi0->Branch("shr2_dot",&_shr2_dot,"shr2_dot/F");
    }

   _p0 = 0;
   _p1 = 0; 
   _p2 = 0;
   _tot_clus = 0;
   _event = -1;

   _purity_v.resize(3) ;
   _complete_v.resize(3) ;

   _cw_purity_v.resize(3) ;
   _cw_complete_v.resize(3) ;

    return true;
  }
  
  bool ShowerCompleteness::analyze(storage_manager* storage) {

    _event++;
    //std::cout<<"EVENT "<<std::endl ;

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

   // Keep track of the charge-weighted hit count as well  
   std::map<int,float> tot_mc_cw_hits_v ; 

   _mc_hit_map.clear();

   // Fill map with hits from mccluster : clusterID
   for( int i = 0; i < ass_mcclus_v.size(); i ++ ){

     //std::cout<<"MCCluster id + hits : "<<i<<", "<<ass_mcclus_v[i].size()<<std::endl ;
     for ( int j = 0; j < ass_mcclus_v[i].size(); j++ ){

       auto hid = ass_mcclus_v[i][j];
       _mc_hit_map[hid] = i ; 

       auto h = ev_hit->at(hid);

       if ( tot_mc_cw_hits_v.find(i) == tot_mc_cw_hits_v.end() )
         tot_mc_cw_hits_v[i] = h.Integral() ;
       else
         tot_mc_cw_hits_v[i] += h.Integral() ;
        
     }
   }

   // Now we have the mc hit info in a map, compute purity and completeness

   // Get the association from cluster -> hit
   auto const & ev_clus = storage->get_data<event_cluster>("ImageClusterHit");
   auto const & ev_ass_c = storage->get_data<larlite::event_ass>("ImageClusterHit");
   auto const & ass_imageclus_v = ev_ass_c->association(ev_clus->id(), ev_hit->id());

   // Get the association from shower -> cluster
   auto const & ev_s = storage->get_data<event_shower>("pi0_candidate_showers"); 
   auto ev_ass_s = storage->get_data<larlite::event_ass>("pi0_candidate_showers");
   auto const& ass_showerreco_v = ev_ass_s->association(ev_s->id(), ev_clus->id());

   // Get mcshower info
   auto const & ev_mcshr = storage->get_data<event_mcshower>("mcreco");
   auto ev_mctruth= storage->get_data<event_mctruth>("generator");

   if ( !ev_clus || ev_clus->size() == 0 || !ev_ass_c || ev_ass_c->size() == 0 ){
     std::cout<<"No cluster..." <<std::endl; 
     return false;
   }   

   if( !ev_s || !ev_s->size() ) { 
     std::cout<<"No shower..." <<std::endl;
     return false;
   }   

   if( !ev_mcshr || !ev_mcshr->size() ) { 
     std::cout<<"No mcshowers..." <<std::endl;
     return false;
   }   

   if(!ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Event has no mctruth info "<<std::endl; 
      return false;
   }

    // First get truth information so we can select only true signal from this final sample
    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    auto traj = nu.Nu().Trajectory();
    auto xvtx = traj.at(traj.size() - 1).X();
    auto yvtx = traj.at(traj.size() - 1).Y();
    auto zvtx = traj.at(traj.size() - 1).Z();

    // Now deal with shower comparisons
   std::vector<int> shr_ids;

   for ( int si = 0; si < ev_mcshr->size(); si++){

     auto s = ev_mcshr->at(si);
     auto st = s.Start();
     auto dist = sqrt( pow(st.X() - xvtx,2) + pow(st.Y() - yvtx,2) + pow(st.Z() - zvtx,2) );

     if ( dist < 0.01 && s.DetProfile().E() > 0 && s.MotherPdgCode() == 111 ){ shr_ids.emplace_back(si) ; }
   }

   if (shr_ids.size() != 2 ) { 
     std::cout<<"More than 2 pi0 showers / pi0s; exiting..."<<std::endl; 
     return false;
   }

   auto s1 = ev_s->at(0) ;
   auto s2 = ev_s->at(1) ;

    std::vector<std::pair<int,int>> mc_to_reco ;
    float max_dot = -1e9;
    int max_mcs(-1), max_recos(-1);

    // Match showers
    for( auto const & mc_id : shr_ids ){
      auto mcs_i = ev_mcshr->at(mc_id);
      auto mag_mcs = sqrt( pow(mcs_i.DetProfile().Px(),2) + pow(mcs_i.DetProfile().Py(),2) + pow(mcs_i.DetProfile().Pz(),2) );

      for( int reco_id = 0; reco_id < 2; reco_id ++ ){

        auto recos_i = ev_s->at(reco_id) ;

        auto mag_reco = sqrt( pow(recos_i.Direction().Px(),2) + pow(recos_i.Direction().Py(),2) + pow(recos_i.Direction().Pz(),2) );
        auto dot = mcs_i.DetProfile().Px() * recos_i.Direction().Px() +
                   mcs_i.DetProfile().Py() * recos_i.Direction().Py() +
                   mcs_i.DetProfile().Pz() * recos_i.Direction().Pz() ;
        dot /= ( mag_mcs * mag_reco );
        if ( fabs(dot) > 1 ) 
          std::cout<<"DOT ! " <<dot <<std::endl ;

        if ( dot > max_dot){
          max_dot = dot;
          max_mcs = mc_id;
          max_recos = reco_id;
         }   
      }   
    }   

    _shr1_dot = max_dot ;
    
    int min_mcs = ( max_mcs == shr_ids[0] ? shr_ids[1] : shr_ids[0] ) ; 
    int min_recos = ( max_recos == 0 ? 1 : 0 ) ; 

    auto mcs2 = ev_mcshr->at(min_mcs) ;
    auto recos2 = ev_s->at(min_recos) ;

    auto mag_mcs = sqrt( pow(mcs2.DetProfile().Px(),2) + pow(mcs2.DetProfile().Py(),2) + pow(mcs2.DetProfile().Pz(),2) );
    auto mag_reco = sqrt( pow(recos2.Direction().Px(),2) + pow(recos2.Direction().Py(),2) + pow(recos2.Direction().Pz(),2) );
    _shr2_dot = mcs2.DetProfile().Px() * recos2.Direction().Px() +
                    mcs2.DetProfile().Py() * recos2.Direction().Py() +
                    mcs2.DetProfile().Pz() * recos2.Direction().Pz() ;
    _shr2_dot /= ( mag_mcs * mag_reco );

    //std::cout<<"SHR 1 + 2 : "<<shr1_dot<<", "<<shr2_dot<<std::endl ;

   // for each reco cluster, find the origin of all hits and calc purity/completeness 
   std::vector<int> pur_ctr_v ;
   std::vector<float> cw_pur_ctr_v ;

   _mass = -10;
   _min_pur = 1e12; 
   _min_com = 1e12; 

   // Loop over showers
   for (size_t i = 0; i < ass_showerreco_v.size(); i++ ){
     
     _tot_clus += ass_showerreco_v.at(i).size() ;

     // Loop over clusters associated to this shower
     for (size_t j = 0; j < ass_showerreco_v.at(i).size(); j++ ){

         auto clus_id = ass_showerreco_v.at(i).at(j); 
         auto iclus = ev_clus->at(clus_id);
     
         int plane = iclus.Plane().Plane ;

         pur_ctr_v.clear();
         cw_pur_ctr_v.clear();

         pur_ctr_v.resize(ass_mcclus_v.size(),0) ;
         cw_pur_ctr_v.resize(ass_mcclus_v.size(),0) ;

         int max_hits = -1;
         int max_cw_hits = -1;
         int max_cid = -1 ;
         float tot_reco_cw_hits = 0;
         
         // Loop through all hits associared to the cluster 
         for ( int k = 0; k < ass_imageclus_v.at(clus_id).size(); k++ ){

           auto hid = ass_imageclus_v.at(clus_id).at(k) ; 
           auto h = ev_hit->at(hid);
           tot_reco_cw_hits += h.Integral() ;
           
           if ( _mc_hit_map.find(hid) != _mc_hit_map.end() ){

             auto mcclus_id = _mc_hit_map[hid] ;

             pur_ctr_v[mcclus_id]++ ; 
             cw_pur_ctr_v[mcclus_id] += h.Integral() ; 

             if( pur_ctr_v[ mcclus_id] > max_hits ){
              
	       //std::cout<<"Reco cluster id + hits: "<<max_cid<<", "<<max_hits<<std::endl;  
               max_hits = pur_ctr_v[mcclus_id];
               max_cid = mcclus_id ; 
               max_cw_hits = cw_pur_ctr_v[mcclus_id] ;
             }
           }
         }

	 float purity = 0. ;
	 float complete = 0. ;

	 float cw_purity = 0. ;
	 float cw_complete = 0.;


         if ( max_cid != -1 ){

           auto tot_mc_hits =  ass_mcclus_v[max_cid].size(); 
           auto tot_reco_hits = ass_imageclus_v[clus_id].size();
           
           purity   = float(max_hits) / tot_reco_hits ;
           complete = float(max_hits) / tot_mc_hits ;

           cw_purity   = float(max_cw_hits) / tot_reco_cw_hits ;
           cw_complete = float(max_cw_hits) / tot_mc_cw_hits_v[max_cid]; 

           if ( (cw_purity < 0.5 || cw_complete < 0.5) && ( _shr1_dot < 0.97 || _shr2_dot < 0.97 ) ){
	     std::cout<<"Event:  "<<_event<<", Plane: "<<plane<<", Pur: "<<cw_purity<<", Com:"<<cw_complete
                      <<", shr1: "<<_shr1_dot<<", "<<_shr2_dot<<std::endl ;

            if( std::find(shit_v.begin(),shit_v.end(),_event) == shit_v.end() )
              shit_v.emplace_back(_event) ; 
           }

	   if ( purity > 1 ){ 
             std::cout<<"*****************************"<<std::endl;
	     std::cout<<"Event: "<<_event<<", clusid: "<<max_cid<<", plane: "<<plane<<std::endl ;
	     std::cout<<"Purtiy : "<<purity<<", "<<complete<<std::endl ;
	     std::cout<<"Max hits, Reco Tot, MC Tot : "<<max_hits<<", "<<tot_reco_hits<<", "<<tot_mc_hits<<std::endl ;
           }
	}

	else{ 
	  std::cout<<"Event:  "<<_event<<", Plane: "<<plane<<", Pur: "<<cw_purity<<", Com:"<<cw_complete<<std::endl ;
          if( std::find(shit_v.begin(),shit_v.end(),_event) == shit_v.end() )
            shit_v.emplace_back(_event) ; 
        }
       

         if( purity < _min_pur )
           _min_pur = purity ;

         if( complete < _min_com )
           _min_com = complete ;

         _purity_v[plane].emplace_back(purity) ;
         _complete_v[plane].emplace_back(complete) ;

         _cw_purity_v[plane].emplace_back(cw_purity) ;
         _cw_complete_v[plane].emplace_back(cw_complete) ;

	 if( plane == 0 )
	   _p0++ ;
	 else if ( plane == 1) 
	   _p1++ ;
	 else if ( plane == 2)
	   _p2++ ;

	//_tot_clus ++ ;
      }
   }

   auto const& shr1 = ev_s->at(0);
   auto const& shr2 = ev_s->at(1);

   // Calc the Opening angle of the showers
   double oangle = acos( shr1.Direction().Dot(shr2.Direction())) ;
   _mass      = sqrt(2 * shr1.Energy() * shr2.Energy() *(1.-cos(oangle)));

   _pi0->Fill(); 

    return true;
  }

  bool ShowerCompleteness::finalize() {

    _tree->Fill();

    std::cout<<" Total fills: "<<_tot_clus<<std::endl ;
    std::cout<<"By plane : "<<_p0<<", "<<_p1<<", "<<_p2<<std::endl ;
    
    std::cout<<"SHITTY EVENTS : "<<shit_v.size()<<std::endl;

    if( _fout ){
      _fout->cd();
      _tree->Write();
      _pi0->Write();
    }
  
    return true;
  }

}
#endif
