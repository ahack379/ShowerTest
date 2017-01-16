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
      _lin_tree->Branch("is_shower",&_is_shower,"is_shower/B");
      _lin_tree->Branch("length",&_length,"length/F");
    }

    _event = -1;
    _event_list = {  29,39,40,47,48,49,59,62,76,79,100,103,104,106,119,143,149,157,158,160,162,171,184,190,196,208,217,221,225,230,247,268,269,271,286,292,306,315,316,319,320,321,364,371,391,393,397,408,417,419,420,424,429,437,440,441,448,452,460,469,476,484,487,489,493,502,504,506,507,514,517,523,533,542,577,579,583,597,599,612,615,620,627,631,645,655,673,681,689,732,736,750,751,757,762,768,776,786,787,790,800,803,805,812,813,825,827,840,845,846,849,853,854,869,875,882,893,897,905,906,912,919,923,931,932,938,943,946,955,956,964,971,997,1023,1026,1034,1037,1039,1045,1047,1090,1093,1115,1116,1120,1122,1123,1129,1138,1141,1144,1147,1154,1161,1190,1193,1205,1216,1219,1223,1225,1227,1243,1260,1275,1285,1299,1305,1306,1308,1313,1325,1328,1337,1338,1355,1359,1384,1385,1391,1394,1404,1411,1430,1438,1446,1452,1461,1468,1469,1475,1477,1485,1490,1499,1509,1516,1529,1533,1539,1541,1545,1551,1554,1560,1561,1573,1585,1588,1591,1598,1612,1623,1643,1645,1646,1651,1661,1663,1665,1666,1668,1669,1671,1672,1675,1690,1691,1699,1707,1709,1717,1726,1735,1736,1760,1763,1764,1771,1777,1780,1786,1793,1795,1798,1801,1804,1805,1811,1820,1845,1848,1849,1859,1861,1866,1871,1874,1878,1886,1887,1896,1897,1903,1919,1926,1939,1940,1944,1945,1946,1949,1961,1962,1966,1983,1984,1995,2004,2017,2027,2028,2031,2033,2037,2044,2049,2056,2058,2066,2071,2074,2079,2082,2083,2089,2100,2121,2122,2132,2141,2151,2163,2167,2171,2174,2195,2201,2215,2223,2227,2236,2239,2259,2263,2265,2285,2307,2311,2314,2320,2338,2352,2354,2355,2372,2375,2386,2388,2391,2393,2395,2398,2408,2413,2417,2420,2424,2427,2428,2431,2433,2435,2440,2455,2465,2471,2478,2488,2489,2490,2492,2496,2498,2503,2506,2513,2535,2536,2550,2552,2560,2573,2576,2592,2604,2615,2631,2649,2659,2663,2698,2699,2709,2731,2736,2739,2746,2799,2808,2830,2831,2838,2851,2875,2878,2883,2889,2897,2906,2913,2914,2920,2933,2950,2955,2968,2979,2980,2981,2986,2995,2998,2999,3000,3011,3019,3031,3034,3057,3071,3075,3077,3078,3086,3111,3118,3123,3127,3134,3135,3138,3143,3155,3156,3164,3165,3171,3182,3184,3204,3206,3207,3214,3225,3246,3249,3254,3265,3267,3277,3285,3306,3307,3312,3327,3338,3346,3348,3350,3352,3357,3362,3363,3373,3381,3389,3392,3398,3411,3414,3433,3434,3444,3445,3461,3465,3471,3474,3484,3490,3499,3504,3514,3517,3520,3542,3547,3548,3552,3558,3581,3601,3608,3611,3612,3618,3620,3622,3628,3638,3640,3648,3668,3681,3701,3718,3724,3731,3740,3744,3755,3774,3777,3783,3784,3799,3827,3828,3844,3846,3848,3854,3860,3861,3866,3869,3880,3885,3887,3888,3897,3898,3905,3908,3918,3925,3927,3934,3949,3951,3954,3959,3960,3963,3970,3974,3984,3988,4007,4019,4032,4034,4046,4051,4058,4059,4062,4063,4066,4079,4090,4104,4109,4128,4135,4136,4138,4141,4147,4148,4157,4160,4167,4177,4184,4205,4226,4229,4240,4247,4256,4263,4278,4294,4304,4312,4313,4316,4325,4334,4338,4340,4341,4352 };




    return true;
  }

  void SepTrkShrNearVtx::Clear(){
    _lin = -999;
    _tll = -999;
    _nhits = -999;
    _is_shower = false;
    _length = -999;
    }
  
  bool SepTrkShrNearVtx::analyze(storage_manager* storage) {

    _event ++ ;
    std::cout<<"\n\nEvent : "<<_event <<std::endl ;

    if(std::find(_event_list.begin(), _event_list.end(), _event) == _event_list.end())
      return false;

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

    if ( ev_shr->size() < 2 ) return false ;

    auto ev_vtx = storage->get_data<event_vertex>("numuCC_vertex");
    if ( !ev_vtx || !ev_vtx->size() ) return false;

    // Get association to shr => clus 
    auto const& ass_clus_v = ev_ass->association(ev_shr->id(), ev_clus->id());
    if (!ass_clus_v.size()) {
      std::cout << "No ass from shower -> clus! " << std::endl;
      return false;
      }   

    // Get association to clus => hit 
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
    std::map<int,std::vector<int>> shr_to_clus ;

    for( size_t s = 0; s < ev_shr->size(); s++) {
      auto const shr = ev_shr->at(s) ;

      for( size_t i = 0; i < ass_clus_v.at(s).size(); i++){
        auto clus = ev_clus->at(ass_clus_v.at(s).at(i)) ;

        //std::cout<<"Plane: "<<clus.Plane().Plane<<std::endl ;
        auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,clus.Plane().Plane);
        auto vtx_w = vtxWT.w;
        auto vtx_t = vtxWT.t + 800 * geomH->TimeToCm() ;

        auto dist_to_vtx = sqrt( pow( clus.StartWire() * geomH->WireToCm() - vtx_w,2) + 
                                 pow( clus.StartTick() * geomH->TimeToCm() - vtx_t,2) );
	
        //std::cout<<"Cluster info : "<<clus.StartTick()*geomH->TimeToCm()<<", "<<clus.StartWire()*geomH->WireToCm()<<std::endl ;
        // If the distance to the vertex is too big, don't consider this "shower" in our sample
        if ( dist_to_vtx > 3.5 ) continue; 

        clus_ids.emplace_back(ass_clus_v.at(s).at(i));
        plane_ids.emplace_back(clus.Plane().Plane);
        shr_ids.emplace_back(s);

	//shr_to_clus[s].emplace_back(ass_clus_v.at(s).at(i)) ;
	//std::cout<<"Adding cluster index : "<<ass_clus_v.at(s).at(i)<<" for shower: "<<s<<std::endl ;
       } 
       shr_to_clus[s] = clus_ids ;
     }

     std::vector<int> mcshr_ids ;
     std::vector<int> mctrk_ids;
     std::vector<int> used_planes ;
     for ( size_t p = 0; p < plane_ids.size(); p++){ 

       if ( std::find(used_planes.begin(),used_planes.end(), plane_ids[p]) != used_planes.end() ) continue;
       used_planes.emplace_back(plane_ids[p]);
       
       auto vtxWT  = geomH->Point_3Dto2D(vtxXYZ,plane_ids[p]); //it->first);
       auto vtx_w = vtxWT.w;
       auto vtx_t = vtxWT.t;  //MC shower is on same offset as vertex

       //std::cout<<"\nIn Plane : "<<plane_ids[p]<<std::endl ;
       for ( size_t k = 0; k < ev_mcshr->size(); k++){ 

         auto traj = ev_mcshr->at(k).Start();
         auto p_st_2d = geomH->Point_3Dto2D(traj.X(),traj.Y(),traj.Z(),plane_ids[p]); //it->first); 
	 auto dist_to_vtx = sqrt( pow( p_st_2d.w - vtx_w,2) + 
                                  pow( p_st_2d.t - vtx_t,2) );
	
         if( dist_to_vtx > 3.5 ) continue; 
	 //std::cout<<"Dist : "<<dist_to_vtx <<std::endl ;

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

     //std::cout<<"Size of reco : "<<shr_ids.size()<<", mc shower: "<<mcshr_ids.size()<<", mc_track: "<<mctrk_ids.size()<<std::endl; 

     // Create map to hold largest scores; maximizing dot product to match mc to reco shower
     std::multimap<float,std::pair<int,int>> score_to_match ;
     std::vector<bool> isShower;
     int prev_id = -1 ;

     //std::cout<<"Shower id sieZ: "<<shr_ids.size()<<std::endl ;

     for( int ii = 0 ; ii < shr_ids.size(); ii++ ){

       //Store shower IDs with cluster pair projections; want to only consider shower once
       if( shr_ids[ii] == prev_id ) continue; 
     
       auto shr_ii = ev_shr->at(shr_ids[ii]);
       auto reco_st = shr_ii.Direction();
       auto mag_ii = sqrt(pow(reco_st.Px(),2) + pow(reco_st.Py(),2) + pow(reco_st.Pz(),2) );
       prev_id = shr_ids[ii];

       for( int jj = 0 ; jj < mcshr_ids.size()+mctrk_ids.size(); jj++ ){
         // Store shower ids for each cluster in each plane; so will have shower duplicates
         if( jj < mcshr_ids.size() ){

           auto mcshr_jj = ev_mcshr->at(mcshr_ids[jj]); 
	   auto mc_st = mcshr_jj.Start() ;

	   auto mag_jj = sqrt(pow(mc_st.Px(),2) + pow(mc_st.Py(),2) + pow(mc_st.Pz(),2) );
	   auto dot = (mc_st.Px()*reco_st.Px() + mc_st.Py()*reco_st.Py() + mc_st.Pz() * reco_st.Pz())/mag_jj/mag_ii ;

	   //std::cout<<"SHOWER DOt is: "<<dot<<", i and j : "<<ii<<", "<<jj<<std::endl ;
	   std::pair<int,int> recomc = { ii, jj }; 
	   score_to_match.emplace(1./dot,recomc) ;//(std::make_pair<ii,jj>));
	   isShower.push_back(true);
	 }
         else{
	   auto jj_adj = jj - mcshr_ids.size() ;
           auto mctrk_jj = ev_mctrk->at(mctrk_ids[jj_adj]); 
	   auto mc_st = mctrk_jj.Start() ;

	   auto mag_jj = sqrt(pow(mc_st.Px(),2) + pow(mc_st.Py(),2) + pow(mc_st.Pz(),2) );
	   auto dot = (mc_st.Px()*reco_st.Px() + mc_st.Py()*reco_st.Py() + mc_st.Pz() * reco_st.Pz())/mag_jj/mag_ii ;

	   //std::cout<<"TRACK DOT is: "<<dot<<std::endl ;
	   std::pair<int,int> recomc = { shr_ids[ii], jj }; 
	   score_to_match.emplace(1./dot,recomc) ;//(std::make_pair<ii,jj>));
	   //final_planes.emplace_back(plane_ids[ii]);
	   isShower.push_back(false);
           } 
         }
       }

       std::vector<int> used_reco_ids ;
       std::vector<int> used_mc_ids ;
       std::vector<std::pair<int,bool>> reco_of_interest ;

       int t = 0;
       for( auto const & i : score_to_match){
         
         if ( std::find(used_reco_ids.begin(),used_reco_ids.end(), i.second.first) != used_reco_ids.end() 
           || std::find(used_mc_ids.begin(),used_mc_ids.end(), i.second.second) != used_mc_ids.end() ) continue;

         //std::cout<<"Score: "<<i.first<<", reco + mc id: "<<i.second.first<<", "<<i.second.second<<std::endl ;

	 used_reco_ids.emplace_back(i.second.first);
	 used_mc_ids.emplace_back(i.second.second);

         auto temp = std::make_pair(i.second.first,isShower[t]);
	 reco_of_interest.emplace_back(temp);
	 t++; 

	 }
    //std::cout<<"Number of matches ! "<<reco_of_interest.size()<<std::endl ;

    for ( size_t j = 0; j < reco_of_interest.size(); j++ ){
      
      auto reco_index = reco_of_interest[j].first ;
      //std::cout<<"Reco shower index : "<<reco_index<<", "<<shr_to_clus[reco_index].size() <<std::endl ;

      for ( size_t k = 0; k < shr_to_clus[reco_index].size(); k++ ){ 

        Clear();

        auto clus_index = shr_to_clus[reco_index].at(k) ;
	//std::cout<<"CLUS INDEX : "<<clus_index <<" for reco shower: "<<reco_index<<std::endl ;
	auto hit_v = ass_hit_v.at(clus_index) ;

        std::vector<double> wire_v ; 
        std::vector<double> time_v ; 

        _nhits = hit_v.size(); 

        for( size_t i = 0; i < hit_v.size(); i++){
	  auto hit_index = hit_v.at(i) ;
          auto h = ev_hit->at(hit_index) ;
 
          wire_v.emplace_back(h.WireID().Wire);
          time_v.emplace_back(h.PeakTime());
	  }

	auto clus = ev_clus->at(clus_index) ;
        _length = sqrt( pow( (clus.EndTick() - clus.StartTick())*geomH->TimeToCm(),2) + 
	                pow( (clus.EndWire() - clus.StartWire())*geomH->WireToCm(),2));

        ::Linearity linear(wire_v,time_v);
        linear._r = 10 ;
   
        linear.linearity(wire_v,time_v) ;
        linear.local_linearity(wire_v,time_v) ;
    
        _lin = linear._lin ;
        _tll = linear._local_lin_truncated_avg;

        for ( auto const & p : reco_of_interest ){
          //std::cout<<"REco of interest "<<p.first<<", "<<j<<std::endl ;
          if ( p.first == j ){
            _is_shower = p.second ; 
            //std::cout<<"j: "<<j<<", Is shower? "<<_is_shower<<std::endl ;
            break ;
            }
          }

        _lin_tree->Fill();
        }
      }

    return true;
  }

  bool SepTrkShrNearVtx::finalize() {

    if(_fout) { _fout->cd(); _lin_tree->Write(); }
  
    return true;
  }

}
#endif
