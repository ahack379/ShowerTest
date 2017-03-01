#ifndef LARLITE_MCVARIABLESTUDY_CXX
#define LARLITE_MCVARIABLESTUDY_CXX

#include "MCVariableStudy.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include <algorithm>

namespace larlite {

  bool MCVariableStudy::initialize() {
    
    if( !_pi0_tree ){
      _pi0_tree = new TTree("pi0_tree","");
      _pi0_tree->Branch("true_pi0_e",&_true_pi0_e,"true_pi0_e/F");
      _pi0_tree->Branch("true_angle",&_true_angle,"true_angle/F");
      _pi0_tree->Branch("true_asym",&_true_asym,"true_asym/F");
      _pi0_tree->Branch("reco_pi0_e",&_reco_pi0_e,"reco_pi0_e/F");
      _pi0_tree->Branch("true_pi0_mom",&_true_pi0_mom,"true_pi0_mom/F");
      _pi0_tree->Branch("true_gamma_e_min",&_true_gamma_e_min,"true_gamma_e_min/F");
      _pi0_tree->Branch("true_gamma_e_max",&_true_gamma_e_max,"true_gamma_e_max/F");
      _pi0_tree->Branch("true_RL_maxE",&_true_RL_maxE,"true_RL_maxE/F");
      _pi0_tree->Branch("true_RL_minE",&_true_RL_minE,"true_RL_minE/F");
      _pi0_tree->Branch("true_nu_e",&_true_nu_e,"true_nu_e/F");
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

  // Cosmics + BNB (With Mesons)
 _event_list = { 29, 39, 48, 59, 62, 76, 79, 100, 103, 106, 119, 157, 160, 171, 184, 190, 196, 208, 217, 269, 286, 292, 306, 316, 319, 321, 364, 371, 391, 397, 419, 424, 429, 437, 440, 448, 452, 460, 484, 487, 489, 493, 502, 504, 506, 514, 517, 533, 542, 577, 597, 612, 615, 620, 655, 673, 681, 732, 750, 751, 762, 768, 776, 786, 790, 800, 803, 840, 846, 849, 853, 869, 882, 893, 897, 912, 919, 931, 943, 964, 971, 1034, 1037, 1045, 1047, 1090, 1093, 1120, 1122, 1123, 1129, 1138, 1141, 1144, 1147, 1161, 1190, 1193, 1216, 1219, 1225, 1243, 1260, 1275, 1285, 1308, 1313, 1325, 1328, 1359, 1384, 1385, 1391, 1394, 1404, 1411, 1430, 1438, 1446, 1452, 1468, 1475, 1485, 1490, 1509, 1516, 1533, 1541, 1545, 1554, 1561, 1573, 1591, 1612, 1643, 1645, 1646, 1651, 1661, 1666, 1671, 1675, 1691, 1707, 1709, 1726, 1735, 1736, 1763, 1771, 1780, 1793, 1801, 1805, 1820, 1845, 1859, 1861, 1886, 1887, 1896, 1897, 1919, 1940, 1944, 1945, 1949, 1962, 1966, 1983, 1984, 1995, 2027, 2028, 2037, 2044, 2049, 2056, 2058, 2071, 2082, 2121, 2122, 2132, 2141, 2151, 2163, 2171, 2174, 2195, 2215, 2223, 2236, 2239, 2259, 2263, 2265, 2285, 2307, 2311, 2320, 2338, 2352, 2354, 2355, 2372, 2386, 2391, 2395, 2398, 2413, 2417, 2420, 2424, 2428, 2435, 2465, 2471, 2478, 2489, 2496, 2498, 2506, 2536, 2560, 2573, 2576, 2604, 2615, 2659, 2698, 2709, 2731, 2736, 2739, 2799, 2808, 2830, 2875, 2883, 2897, 2913, 2920, 2933, 2980, 2981, 2986, 2995, 3000, 3019, 3031, 3057, 3071, 3075, 3077, 3078, 3118, 3123, 3127, 3134, 3135, 3143, 3155, 3164, 3165, 3171, 3204, 3206, 3207, 3225, 3246, 3249, 3265, 3267, 3306, 3312, 3327, 3346, 3348, 3350, 3362, 3381, 3392, 3411, 3414, 3433, 3434, 3484, 3490, 3504, 3517, 3520, 3542, 3547, 3552, 3558, 3581, 3601, 3611, 3612, 3622, 3638, 3668, 3718, 3740, 3744, 3783, 3784, 3827, 3828, 3844, 3846, 3848, 3860, 3861, 3866, 3880, 3887, 3888, 3898, 3905, 3908, 3918, 3927, 3949, 3951, 3954, 3959, 3963, 3974, 3984, 4034, 4046, 4062, 4066, 4090, 4104, 4135, 4136, 4138, 4141, 4148, 4157, 4184, 4205, 4226, 4229, 4240, 4256, 4263, 4312, 4313, 4334, 4341, 4352 };

 // BNB Only (No Mesons)
 //_event_list = { 0, 10, 13, 16, 21, 23, 87, 102, 108, 135, 146, 147, 149, 150, 166, 173, 184, 197, 215, 250, 264, 276, 291, 316, 320, 366, 372, 373, 379, 381, 384, 428, 434, 437, 449, 450, 490, 555, 574, 611, 613, 624, 635, 638, 653, 654, 664, 666, 687, 694, 714, 746, 756, 757, 758, 795, 809, 834, 838, 848, 849, 881, 886, 902, 920, 972, 1001, 1060, 1076, 1077, 1112, 1119, 1139, 1160, 1162, 1212, 1216, 1219, 1220, 1248, 1271, 1286, 1297, 1321, 1335, 1339, 1343, 1376, 1395, 1409, 1412, 1423, 1430, 1434, 1468, 1480, 1528, 1609, 1647, 1657, 1661, 1691, 1712, 1725, 1747, 1807, 1812, 1837, 1846, 1852, 1856, 1860, 1877, 1888, 1899, 1987, 1988, 1996, 1998, 2011, 2060, 2068, 2099, 2120, 2128, 2143, 2152, 2215, 2242, 2255, 2259, 2271, 2291, 2298, 2328, 2351, 2374, 2383, 2405, 2419, 2421, 2425, 2446, 2462, 2495, 2502, 2515, 2527, 2557, 2573, 2580, 2608, 2647, 2652, 2655, 2664, 2667, 2681, 2682, 2692, 2716, 2718, 2719, 2731, 2735, 2746, 2757, 2762, 2774, 2795, 2807, 2822, 2836, 2840, 2853, 2857, 2864, 2868, 2877, 2889, 2907, 2908, 2959, 2971, 2982, 2995, 3011, 3028, 3030, 3060, 3069, 3079, 3111, 3130, 3139, 3146, 3157, 3167, 3171, 3178, 3210, 3216, 3223, 3224, 3245, 3266, 3275, 3287, 3297, 3318, 3336, 3370, 3401, 3409, 3438, 3465, 3485, 3495, 3525, 3530, 3589, 3619, 3641, 3660, 3667, 3670, 3675, 3677, 3685, 3687, 3719, 3720, 3752, 3754, 3780, 3789, 3815, 3816, 3825, 3826, 3877, 3898, 3912, 3935, 3963, 4024, 4034, 4057, 4129, 4130, 4135, 4173, 4189, 4191, 4208, 4272, 4275, 4298, 4301, 4367, 4372, 4375, 4389, 4402, 4403, 4429, 4435, 4448, 4449, 4464, 4500, 4505, 4527, 4568, 4580, 4588, 4604, 4611, 4614, 4617, 4639, 4655, 4664, 4668, 4671, 4695, 4702, 4712, 4723, 4735, 4786, 4792, 4821, 4825, 4835, 4846, 4848, 4849, 4864, 4891, 4898, 4957, 5005, 5029, 5051, 5069, 5093, 5101, 5125, 5158, 5177, 5199, 5221, 5238, 5284, 5329, 5392, 5399, 5400, 5402, 5429, 5437, 5457, 5471, 5475, 5489, 5491, 5518, 5578, 5606, 5621, 5667, 5679, 5683, 5696, 5716, 5717, 5718, 5725, 5760, 5762, 5778, 5814, 5843, 5849, 5854, 5883, 5895, 5919, 5934, 5968, 6009, 6059 };
  
   std::cout<<"Event list size: "<<_event_list.size() <<std::endl ;

    return true;
  }

  void MCVariableStudy::Clear(){
    _true_pi0_e = -999;
    _true_angle = -9;
    _true_asym = -9;
    _reco_pi0_e = -999;
    _true_gamma_e_max = 0 ;
    _true_gamma_e_min = 0 ;
    _true_pi0_mom = -999;
    _true_nu_e = -999;
    _true_RL_maxE = 0;
    _true_RL_minE = 0;

    _n_true_pi0 = 0;

    _true_gamma_e = -9;
    _reco_gamma_e = -9;
    _true_rad_l = -9;
    _reco_rad_l = -9;
    _true_reco_dot = -9;
    }
  
  bool MCVariableStudy::analyze(storage_manager* storage) {

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

    auto & nu  = ev_mctruth->at(0).GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    if( xyz[0] < 20 || xyz[0] > 236.25 || xyz[1] < -96.5 || xyz[1] > 96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
      return false ;

    auto parts = ev_mctruth->at(0).GetParticles();

    int n_mes = 0;
    int n_lep = 0;
    int n_mu = 0;
    
    std::vector<float> start(3,0) ;
    
    for ( auto const & p : parts ){
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){
        _n_true_pi0 ++;
        _true_pi0_e = p.Trajectory().at(0).E()*1000; 
        start[0] = p.Trajectory().at(0).X();
        start[1] = p.Trajectory().at(0).Y();
        start[2] = p.Trajectory().at(0).Z();
	_true_pi0_mom = sqrt( pow(p.Trajectory().at(0).Px(),2) + pow(p.Trajectory().at(0).Py(),2) + pow(p.Trajectory().at(0).Pz(),2) )*1000; 
        }

      if( p.StatusCode() == 1 && p.PdgCode() == 13 )
        n_mu ++; 

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 211 || abs(p.PdgCode()) == 321 ||  p.PdgCode() == 130 
                                   || p.PdgCode() == 310 || abs(p.PdgCode()) == 311 ) ){
        n_mes ++; 
        }   

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 11 || p.PdgCode() == -13))
        n_lep ++; 

      }

    if( n_mu != 1 || _n_true_pi0 != 1) return false ; // || n_lep != 0 || n_mes != 0) return false;
    //if ( _n_true_pi0 != 1 ) return false; 

    _true_nu_e = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().at(0).E() ;
    std::cout<<"TRue nu energy: "<<_true_nu_e<<std::endl; 
      
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

      auto d1 = ev_mcs->at(shr_ids[0]).DetProfile();
      auto d2 = ev_mcs->at(shr_ids[1]).DetProfile();

      auto mag1 = sqrt( s1.Px()*s1.Px()+s1.Py()*s1.Py()+s1.Pz()*s1.Pz() );
      auto mag2 = sqrt( s2.Px()*s2.Px()+s2.Py()*s2.Py()+s2.Pz()*s2.Pz() );
      auto dot = s1.Px()*s2.Px() + s1.Py()*s2.Py() + s1.Pz()*s2.Pz() ;

      auto radL1 = sqrt( pow(d1.X() - xyz[0],2) + pow(d1.Y() - xyz[1],2) + pow(d1.Z() - xyz[2],2) ); 
      auto radL2 = sqrt( pow(d2.X() - xyz[0],2) + pow(d2.Y() - xyz[1],2) + pow(d2.Z() - xyz[2],2) ); 
      std::cout<<"WOAH :" <<s1.X()<<", "<<xyz[0]<<std::endl;

      float e1 = s1.E() ;
      float e2 = s2.E() ;

      //std::cout<<"DOT: "<<dot<<", "<<mag1<<", "<<mag2<<std::endl ; 

      _true_gamma_e_min = e1 < e2 ? e1 : e2 ;
      _true_gamma_e_max = e1 > e2 ? e1 : e2 ;

      _true_RL_minE = e1 < e2 ? radL1 : radL2 ;
      _true_RL_maxE = e1 > e2 ? radL1 : radL2 ;

      _true_asym = e1 > e2 ? e2/e1 : e1/e2 ;
      _true_angle = acos( dot / mag1 / mag2 ); 
    }

   _pi0_tree->Fill();
      
   // auto ev_s = storage->get_data<event_shower>("showerreco"); 
   // if( !ev_s || !ev_s->size() ){
   //   std::cout<<"No shower..." <<std::endl;
   //   return false;
   //   }

   // if( ev_s->size() < 2 ) return false ;
   // 
   // float min_IP = 1e9;
   // int min_it = -1 ;

   // std::vector<std::pair<int,int>> CandidatePairs;
   // std::vector<int> cand_ids;

    //std::cout<<ev_s->size()<<" Showers : "<<std::endl; 

   // for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

   //         auto const& shr1 = ev_s->at(s1);

   //         for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

   //             if (s2 <= s1) continue;

   //             auto const& shr2 = ev_s->at(s2);

   //             geoalgo::Vector_t dflipA(-1.*shr1.Direction()) ; 
   //             geoalgo::Vector_t dflipB(-1.*shr2.Direction()) ;

   //             // Make the backwards projection for the showers
   //             auto bksa = ::geoalgo::HalfLine_t(shr1.ShowerStart(),dflipA);
   //             auto bksb = ::geoalgo::HalfLine_t(shr2.ShowerStart(),dflipB);

   //             // Calc the Opening angle of the showers
   //             double Oangle = acos( shr1.Direction().Dot(shr2.Direction())) ; 

   //             // Calc the vertex point of the two showers. the true designated backwards project
   //             geoalgo::Point_t vertex(3);
   //             
   //             auto st1 = shr1.ShowerStart();
   //             auto st2 = shr2.ShowerStart();
   //             auto dir1 = shr1.Direction();
   //             auto dir2 = shr2.Direction();
   //             geoalgo::HalfLine_t shr1_hl(st1.X(),st1.Y(),st1.Z(),dir1.X(),dir1.Y(), dir1.Z() );
   //             geoalgo::HalfLine_t shr2_hl(st2.X(),st2.Y(),st2.Z(),dir2.X(),dir2.Y(), dir2.Z() );

   //             _geoAlgo.commonOrigin(shr1_hl, shr2_hl, vertex, true);

   //             // Calc Diretion of two correlated shower
   //             geoalgo::Vector_t momentum(3);// need to fill out
   //             geoalgo::Vector_t mom_vect(shr2.Direction()*shr1.Energy(2) +shr1.Direction()*shr2.Energy(2)) ;
   //             mom_vect.Normalize();
   //             momentum = mom_vect * sqrt(pow(shr1.Energy(),2)+pow(shr2.Energy(),2)+2*shr2.Energy()*shr1.Energy()*cos(Oangle));

   //             //===========================================
   //             auto _RIP = pow(_geoAlgo.SqDist(bksa,bksb),0.5);
   //             auto _RadL_A = vertex.Dist(shr1.ShowerStart());
   //             auto _RadL_B = vertex.Dist(shr2.ShowerStart());

   //             if(Oangle < 0.35){
   //               //std::cout<<"Bad Angle: "<< Oangle<<std::endl ; 
   //               continue;
   //              }

   //             if(pow(_geoAlgo.SqDist(bksa,bksb),0.5) > 4){
   //                //std::cout<<"Ip? "<<pow(_geoAlgo.SqDist(bksa,bksb),0.5)<<std::endl ; 
   //                continue;
   //               }

   //             if(_RadL_A > 62. || _RadL_B > 62. ){ 
   //                //std::cout<<"Rad Length : "<<_RadL_A<<", "<<_RadL_B<<std::endl ;
   //                continue; 
   //                 }

   //             if( vertex[0] > 256.35 || vertex[0] < 0. || vertex[1] > 116.5  
   //               || vertex[1] < -116.5 || vertex[2] > 1036.8 || vertex[2] < 0. ){

   //               std::cout<<"Vertex out of bounds : "<<vertex[0]<<", "<<vertex[1]<<", "<<vertex[2]<<std::endl ;
   //               continue;
   //               }

   //             // Bunch of cuts
   //             if( _RIP < min_IP ){
   //               min_IP = _RIP ;
   //               min_it = CandidatePairs.size();
   //               }

   //             CandidatePairs.push_back(std::make_pair(s1,s2)); 
   //     	cand_ids.emplace_back(s1);
   //     	cand_ids.emplace_back(s2);

   //         }// shower ID 2 
   //     }// shower ID 1 

   //    std::cout<<"Cand size: "<<cand_ids.size()<<std::endl ;
   //    // Fill trees 
   //    if( cand_ids.size() == 2 && min_it != -1){ 
   // 
   //        auto s1 = ev_s->at(cand_ids[0]) ;
   //        auto s2 = ev_s->at(cand_ids[1]) ;
   //       
   //        _reco_pi0_e = s1.Energy(2) + s2.Energy(2); 
   //        _pi0_tree->Fill();

   //        std::vector<std::pair<int,int>> mc_to_reco ;
   //        float max_dot = -1e9;
   //        int max_mcs(-1), max_recos(-1), min_mcs(-1), min_recos(-1);

   //        // Match showers
   //        for( auto const & mc_id : shr_ids ){
   //          auto mcs_i = ev_mcs->at(mc_id);
   //          auto mag_mcs = sqrt( pow(mcs_i.Start().Px(),2) + pow(mcs_i.Start().Py(),2) + pow(mcs_i.Start().Pz(),2) ); 
   //          
   //          for( auto const & reco_id : cand_ids ){

   //            auto recos_i = ev_s->at(reco_id) ;

   //            auto mag_reco = sqrt( pow(recos_i.Direction().Px(),2) + pow(recos_i.Direction().Py(),2) + pow(recos_i.Direction().Pz(),2) ); 
   //            auto dot = mcs_i.Start().Px() * recos_i.Direction().Px() + 
   //                       mcs_i.Start().Py() * recos_i.Direction().Py() + 
   //     		  mcs_i.Start().Pz() * recos_i.Direction().Pz() ;
   //            dot /= ( mag_mcs * mag_reco );

   //            if ( dot > max_dot){
   //              max_dot = dot;
   //              max_mcs = mc_id;
   //              max_recos = reco_id;
   //             }
   //           } 
   //         }
   //         
   //        min_mcs = ( max_mcs == shr_ids[0] ? shr_ids[1] : shr_ids[0] ) ;
   //        min_recos = ( max_recos == cand_ids[0] ? cand_ids[1] : cand_ids[0] ) ;

   //        auto mcs_1 = ev_mcs->at(max_mcs) ;
   //        auto recos_1 = ev_s->at(max_recos) ;

   //        _true_gamma_e = mcs_1.Start().E();
   //        _reco_gamma_e = recos_1.Energy(2);
   //        _true_reco_dot = max_dot;
   //        _gamma_tree->Fill();
   //        std::cout<<"Reco and tru gamma e: "<<_reco_gamma_e<<", "<<_true_gamma_e<<std::endl ;
   //        std::cout<<"Start of shower: "<<mcs_1.Start().X()<<", "<<mcs_1.Start().Y() <<", "<<mcs_1.Start().Z()<<std::endl ;
   //        std::cout<<"End of shower: "<<mcs_1.End().X()<<", "<<mcs_1.End().Y() <<", "<<mcs_1.End().Z()<<std::endl ;

   //        auto mcs_2 = ev_mcs->at(min_mcs) ;
   //        auto recos_2 = ev_s->at(min_recos) ;
   //        auto dot = mcs_2.Start().Px() * recos_2.Direction().Px() + mcs_2.Start().Py() * recos_2.Direction().Py() + mcs_2.Start().Pz() * recos_2.Direction().Pz() ;
   //        auto mag_reco = sqrt( pow(recos_2.Direction().Px(),2) + pow(recos_2.Direction().Py(),2) + pow(recos_2.Direction().Pz(),2) ); 
   //        auto mag_mcs = sqrt( pow(mcs_2.Start().Px(),2) + pow(mcs_2.Start().Py(),2) + pow(mcs_2.Start().Pz(),2) ); 
   //        dot /= (mag_mcs * mag_reco);

   //        _true_gamma_e = mcs_2.Start().E();
   //        _reco_gamma_e = recos_2.Energy(2);
   //        _true_reco_dot = dot;
   //        _gamma_tree->Fill();

   //        std::cout<<"\nReco and tru gamma e: "<<_reco_gamma_e<<", "<<_true_gamma_e<<std::endl ;
   //        std::cout<<"Start of shower: "<<mcs_2.Start().X()<<", "<<mcs_2.Start().Y() <<", "<<mcs_2.Start().Z()<<std::endl ;
   //        std::cout<<"End of shower: "<<mcs_2.End().X()<<", "<<mcs_2.End().Y() <<", "<<mcs_2.End().Z()<<std::endl ;
   //      }
  
    return true;
  }

  bool MCVariableStudy::finalize() {

    if(_fout) { 
      _fout->cd(); 
      _pi0_tree->Write(); 
      _gamma_tree->Write(); 
      }
  
    return true;
  }

}
#endif
