#ifndef LARLITE_TRUTHDISTRIBS_CXX
#define LARLITE_TRUTHDISTRIBS_CXX

#include "TruthDistribs.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include <algorithm>

namespace larlite {

  bool TruthDistribs::initialize() {
    
    if( !_pi0_tree ){
      _pi0_tree = new TTree("pi0_tree","");
      _pi0_tree->Branch("true_pi0_e",&_true_pi0_e,"true_pi0_e/F");
      _pi0_tree->Branch("true_angle",&_true_angle,"true_angle/F");
      _pi0_tree->Branch("true_asym",&_true_asym,"true_asym/F");
      _pi0_tree->Branch("reco_pi0_e",&_reco_pi0_e,"reco_pi0_e/F");
      _pi0_tree->Branch("event",&_event,"event/I");
      }

    if( !_gamma_tree ){
      _gamma_tree = new TTree("gamma_tree","");
      _gamma_tree->Branch("true_gamma_e",&_true_gamma_e,"true_gamma_pi0_e/F");
      _gamma_tree->Branch("reco_gamma_e",&_reco_gamma_e,"reco_gamma_pi0_e/F");
      _gamma_tree->Branch("true_reco_dot",&_true_reco_dot,"true_reco_dot/F");
      _gamma_tree->Branch("true_rad_l",&_true_rad_l,"true_rad_l/F");
      _gamma_tree->Branch("reco_rad_l",&_reco_rad_l,"reco_rad_l/F");
      }

    _event = -1;
    _event_list = {17,29,30,36,48,49,55,59,62,76,79,102,103,104,108,119,133,150,157,158,159,160,162,171,190,193,196,200,208,217,219,221,225,230,251,268,269,271,283,292,297,306,315,316,319,320,321,331,339,371,391,397,408,419,420,425,429,433,437,441,448,452,460,469,470,483,484,491,500,502,506,514,517,533,542,557,577,579,583,588,597,599,607,612,615,617,618,619,620,627,631,644,645,655,673,681,684,689,709,710,732,736,751,757,762,776,787,790,793,803,805,812,813,825,827,840,846,849,850,853,854,869,875,882,893,896,897,905,907,912,919,921,923,932,938,943,945,946,950,955,956,964,971,976,986,1000,1016,1023,1026,1034,1036,1037,1039,1043,1045,1047,1073,1090,1093,1098,1100,1109,1112,1116,1120,1122,1123,1129,1141,1144,1147,1154,1161,1178,1190,1201,1205,1216,1219,1223,1227,1243,1261,1270,1282,1285,1287,1299,1304,1306,1308,1313,1325,1328,1341,1355,1359,1384,1385,1386,1391,1394,1404,1411,1430,1437,1438,1446,1455,1461,1468,1475,1477,1485,1490,1494,1499,1506,1509,1516,1525,1529,1533,1539,1541,1551,1554,1559,1560,1561,1573,1585,1588,1591,1612,1623,1625,1630,1643,1645,1661,1665,1666,1669,1671,1672,1675,1690,1699,1704,1707,1709,1717,1726,1735,1736,1754,1760,1763,1764,1771,1777,1780,1786,1793,1795,1804,1805,1809,1811,1813,1820,1825,1845,1848,1859,1866,1871,1874,1877,1878,1887,1896,1897,1906,1919,1926,1939,1940,1944,1945,1949,1961,1962,1966,1968,1974,1984,1990,1991,1995,2017,2021,2022,2027,2028,2033,2044,2049,2056,2058,2062,2066,2067,2074,2081,2082,2083,2089,2093,2096,2100,2116,2122,2132,2139,2141,2151,2157,2163,2167,2174,2195,2199,2201,2212,2213,2215,2219,2236,2237,2239,2259,2263,2265,2276,2285,2290,2300,2307,2311,2314,2320,2338,2343,2354,2355,2372,2386,2388,2393,2398,2413,2417,2418,2420,2424,2427,2428,2431,2435,2438,2465,2468,2471,2478,2481,2489,2490,2492,2498,2503,2506,2513,2535,2536,2541,2550,2552,2554,2560,2573,2576,2580,2592,2604,2607,2615,2625,2647,2659,2666,2673,2698,2699,2701,2707,2710,2720,2731,2736,2739,2766,2778,2792,2799,2810,2811,2819,2825,2831,2834,2842,2873,2875,2878,2883,2906,2914,2919,2920,2927,2933,2949,2950,2955,2968,2980,2981,2985,2986,2995,2999,3000,3011,3019,3022,3030,3031,3034,3051,3052,3057,3071,3075,3088,3103,3111,3113,3116,3118,3123,3127,3133,3135,3138,3143,3154,3155,3156,3164,3165,3167,3171,3178,3182,3204,3206,3207,3225,3249,3250,3254,3255,3256,3267,3285,3293,3296,3303,3306,3307,3312,3335,3338,3346,3348,3357,3362,3373,3389,3394,3398,3413,3414,3434,3444,3446,3461,3465,3471,3474,3483,3484,3490,3499,3514,3520,3538,3542,3547,3548,3552,3558,3578,3581,3586,3595,3601,3604,3611,3612,3613,3618,3628,3638,3640,3663,3681,3691,3704,3718,3724,3755,3771,3774,3784,3799,3817,3827,3828,3831,3844,3846,3860,3861,3862,3863,3866,3869,3880,3884,3887,3888,3898,3902,3905,3906,3908,3916,3925,3927,3934,3937,3949,3951,3959,3963,3970,3974,3980,3984,3988,4000,4007,4011,4029,4032,4034,4042,4044,4051,4058,4059,4064,4079,4098,4100,4101,4104,4109,4128,4135,4136,4138,4141,4144,4147,4148,4157,4174,4177,4184,4205,4226,4229,4232,4240,4243,4247,4278,4294,4304,4312,4325,4340,4341,4352};
  
   std::cout<<"Event list size: "<<_event_list.size() <<std::endl ;

    return true;
  }

  void TruthDistribs::Clear(){
    _true_pi0_e = -999;
    _true_angle = -9;
    _true_asym = -9;
    _reco_pi0_e = -999;

    _n_true_pi0 = 0;

    _true_gamma_e = -9;
    _reco_gamma_e = -9;
    _true_rad_l = -9;
    _reco_rad_l = -9;
    _true_reco_dot = -9;
    }
  
  bool TruthDistribs::analyze(storage_manager* storage) {

    _event ++; 
    Clear() ;

    if(std::find(_event_list.begin(), _event_list.end(), _event) == _event_list.end()) 
      return false;

    std::cout<<"\nEvent is : "<<_event <<std::endl ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ){
      std::cout<<"No mctruth..." <<std::endl;
      return false;
      }

    auto parts = ev_mctruth->at(0).GetParticles();

    bool pi0 = false;
    bool mu  = false ;
    
    std::vector<float> start(3,0) ;
    
    for ( auto const & p : parts ){
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){
        _n_true_pi0 ++;
        _true_pi0_e = p.Trajectory().at(0).E()*1000; 
        start[0] = p.Trajectory().at(0).X();
        start[1] = p.Trajectory().at(0).Y();
        start[2] = p.Trajectory().at(0).Z();
        }
      }   
    
    if ( _n_true_pi0 != 1 ) return false; 
      
   // Replace pi0 energy with combined shower energy 
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco"); 
    if(!ev_mcs || !ev_mcs->size() ){
      std::cout<<"No mcshower..." <<std::endl;
      return false;
      }

    std::vector<int> shr_ids;
     
    for ( int si = 0; si < ev_mcs->size(); si++){ 

      auto s = ev_mcs->at(si);

      if( s.PdgCode() != 22 ) continue; //|| abs(s.MotherPdgCode()) == 13 ) continue;
      
      auto st = s.Start();
      auto dist = sqrt( pow(st.X() - start[0],2) + pow(st.Y() - start[1],2) + pow(st.Z() - start[2],2) );
      
      if ( dist < 0.001 ){
        shr_ids.emplace_back(si) ;
        //std::cout<<"Dist: "<<dist<<", energies shr: "<<ev_mcs->at(si).Start().E()<<std::endl;
	}
    }
    
    if( shr_ids.size() == 2 ){

      auto s1 = ev_mcs->at(shr_ids[0]).Start();
      auto s2 = ev_mcs->at(shr_ids[1]).Start();

      auto mag1 = sqrt( s1.Px()*s1.Px()+s1.Py()*s1.Py()+s1.Pz()*s1.Pz() );
      auto mag2 = sqrt( s2.Px()*s2.Px()+s2.Py()*s2.Py()+s2.Pz()*s2.Pz() );
      auto dot = s1.Px()*s2.Px() + s1.Py()*s2.Py() + s1.Pz()*s2.Pz() ;

      float e1 = s1.E() ;
      float e2 = s2.E() ;

      //std::cout<<"DOT: "<<dot<<", "<<mag1<<", "<<mag2<<std::endl ; 

      _true_asym = e1 > e2 ? e2/e1 : e1/e2 ;
      _true_angle = acos( dot / mag1 / mag2 ); 
    }
      
    auto ev_s = storage->get_data<event_shower>("showerreco"); 
    if( !ev_s || !ev_s->size() ){
      std::cout<<"No shower..." <<std::endl;
      return false;
      }

    if( ev_s->size() < 2 ) return false ;
    
    float min_IP = 1e9;
    int min_it = -1 ;

    std::vector<std::pair<int,int>> CandidatePairs;
    std::vector<int> cand_ids;

    //std::cout<<ev_s->size()<<" Showers : "<<std::endl; 

    for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

            auto const& shr1 = ev_s->at(s1);

            for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

                if (s2 <= s1) continue;

                auto const& shr2 = ev_s->at(s2);

                geoalgo::Vector_t dflipA(-1.*shr1.Direction()) ; 
                geoalgo::Vector_t dflipB(-1.*shr2.Direction()) ;

                // Make the backwards projection for the showers
                auto bksa = ::geoalgo::HalfLine_t(shr1.ShowerStart(),dflipA);
                auto bksb = ::geoalgo::HalfLine_t(shr2.ShowerStart(),dflipB);

                // Calc the Opening angle of the showers
                double Oangle = acos( shr1.Direction().Dot(shr2.Direction())) ; 

                // Calc the vertex point of the two showers. the true designated backwards project
                geoalgo::Point_t vertex(3);
                
                auto st1 = shr1.ShowerStart();
                auto st2 = shr2.ShowerStart();
                auto dir1 = shr1.Direction();
                auto dir2 = shr2.Direction();
                geoalgo::HalfLine_t shr1_hl(st1.X(),st1.Y(),st1.Z(),dir1.X(),dir1.Y(), dir1.Z() );
                geoalgo::HalfLine_t shr2_hl(st2.X(),st2.Y(),st2.Z(),dir2.X(),dir2.Y(), dir2.Z() );

                _geoAlgo.commonOrigin(shr1_hl, shr2_hl, vertex, true);

                // Calc Diretion of two correlated shower
                geoalgo::Vector_t momentum(3);// need to fill out
                geoalgo::Vector_t mom_vect(shr2.Direction()*shr1.Energy(2) +shr1.Direction()*shr2.Energy(2)) ;
                mom_vect.Normalize();
                momentum = mom_vect * sqrt(pow(shr1.Energy(),2)+pow(shr2.Energy(),2)+2*shr2.Energy()*shr1.Energy()*cos(Oangle));

                //===========================================
                auto _RIP = pow(_geoAlgo.SqDist(bksa,bksb),0.5);
                auto _RadL_A = vertex.Dist(shr1.ShowerStart());
                auto _RadL_B = vertex.Dist(shr2.ShowerStart());

                if(Oangle < 0.35){
                  //std::cout<<"Bad Angle: "<< Oangle<<std::endl ; 
                  continue;
                 }

                if(pow(_geoAlgo.SqDist(bksa,bksb),0.5) > 4){
                   //std::cout<<"Ip? "<<pow(_geoAlgo.SqDist(bksa,bksb),0.5)<<std::endl ; 
                   continue;
                  }

                if(_RadL_A > 62. || _RadL_B > 62. ){ 
                   //std::cout<<"Rad Length : "<<_RadL_A<<", "<<_RadL_B<<std::endl ;
                   continue; 
                    }

                if( vertex[0] > 256.35 || vertex[0] < 0. || vertex[1] > 116.5  
                  || vertex[1] < -116.5 || vertex[2] > 1036.8 || vertex[2] < 0. ){

                  std::cout<<"Vertex out of bounds : "<<vertex[0]<<", "<<vertex[1]<<", "<<vertex[2]<<std::endl ;
                  continue;
                  }

                // Bunch of cuts
                if( _RIP < min_IP ){
                  min_IP = _RIP ;
                  min_it = CandidatePairs.size();
                  }

                CandidatePairs.push_back(std::make_pair(s1,s2)); 
		cand_ids.emplace_back(s1);
		cand_ids.emplace_back(s2);

            }// shower ID 2 
        }// shower ID 1 

       std::cout<<"Cand size: "<<cand_ids.size()<<std::endl ;
       // Fill trees 
       if( cand_ids.size() == 2 && min_it != -1){ 
    
           auto s1 = ev_s->at(cand_ids[0]) ;
           auto s2 = ev_s->at(cand_ids[1]) ;
          
           _reco_pi0_e = s1.Energy(2) + s2.Energy(2); 
           _pi0_tree->Fill();

           std::vector<std::pair<int,int>> mc_to_reco ;
	   float max_dot = -1e9;
	   int max_mcs(-1), max_recos(-1), min_mcs(-1), min_recos(-1);

           // Match showers
           for( auto const & mc_id : shr_ids ){
             auto mcs_i = ev_mcs->at(mc_id);
	     auto mag_mcs = sqrt( pow(mcs_i.Start().Px(),2) + pow(mcs_i.Start().Py(),2) + pow(mcs_i.Start().Pz(),2) ); 
	     
	     for( auto const & reco_id : cand_ids ){

	       auto recos_i = ev_s->at(reco_id) ;

	       auto mag_reco = sqrt( pow(recos_i.Direction().Px(),2) + pow(recos_i.Direction().Py(),2) + pow(recos_i.Direction().Pz(),2) ); 
	       auto dot = mcs_i.Start().Px() * recos_i.Direction().Px() + 
	                  mcs_i.Start().Py() * recos_i.Direction().Py() + 
			  mcs_i.Start().Pz() * recos_i.Direction().Pz() ;
	       dot /= ( mag_mcs * mag_reco );

               if ( dot > max_dot){
	         max_dot = dot;
	         max_mcs = mc_id;
	         max_recos = reco_id;
                }
	      } 
	    }
	    
	   min_mcs = ( max_mcs == shr_ids[0] ? shr_ids[1] : shr_ids[0] ) ;
	   min_recos = ( max_recos == cand_ids[0] ? cand_ids[1] : cand_ids[0] ) ;

	   auto mcs_1 = ev_mcs->at(max_mcs) ;
	   auto recos_1 = ev_s->at(max_recos) ;

	   _true_gamma_e = mcs_1.Start().E();
	   _reco_gamma_e = recos_1.Energy(2);
	   _true_reco_dot = max_dot;
	   _gamma_tree->Fill();

	   auto mcs_2 = ev_mcs->at(min_mcs) ;
	   auto recos_2 = ev_s->at(min_recos) ;
	   auto dot = mcs_2.Start().Px() * recos_2.Direction().Px() + mcs_2.Start().Py() * recos_2.Direction().Py() + mcs_2.Start().Pz() * recos_2.Direction().Pz() ;
	   auto mag_reco = sqrt( pow(recos_2.Direction().Px(),2) + pow(recos_2.Direction().Py(),2) + pow(recos_2.Direction().Pz(),2) ); 
	   auto mag_mcs = sqrt( pow(mcs_2.Start().Px(),2) + pow(mcs_2.Start().Py(),2) + pow(mcs_2.Start().Pz(),2) ); 
	   dot /= (mag_mcs * mag_reco);

	   _true_gamma_e = mcs_2.Start().E();
	   _reco_gamma_e = recos_2.Energy(2);
	   _true_reco_dot = dot;
	   _gamma_tree->Fill();
         }
  
    return true;
  }

  bool TruthDistribs::finalize() {

    if(_fout) { 
      _fout->cd(); 
      _pi0_tree->Write(); 
      _gamma_tree->Write(); 
      }
  
    return true;
  }

}
#endif
