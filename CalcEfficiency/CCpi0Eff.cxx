#ifndef LARLITE_CCPI0EFF_CXX
#define LARLITE_CCPI0EFF_CXX

#include "CCpi0Eff.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/vertex.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/LArProperties.h"
#include "LArUtil/TimeService.h"

namespace larlite {

  bool CCpi0Eff::initialize() {    

    _event = -1; 
    _signal = 0;
    _event_list.clear();

    _multpi0 = 0;
    _mesons = 0;
    _leptons = 0;
    _gammas = 0;
    _tot_ccpi0 = 0;
    _nc_pi0 = 0;

    //_diff_list = { 29, 184, 196, 306, 316, 440, 460, 502, 514, 615, 750, 768, 790, 800, 803, 882, 893, 912, 1120, 1122, 1138, 1141, 1193, 1216, 1260, 1313, 1384, 1452, 1468, 1485, 1554, 1591, 1646, 1675, 1691, 1735, 1763, 1771, 1793, 1820, 1887, 1896, 1962, 1966, 2028, 2058, 2071, 2082, 2239, 2307, 2311, 2354, 2355, 2386, 2395, 2465, 2478, 2506, 2536, 2659, 2731, 2736, 2808, 2830, 2883, 2913, 2920, 2980, 3071, 3134, 3155, 3165, 3171, 3348, 3350, 3392, 3411, 3434, 3484, 3490, 3520, 3558, 3622, 3668, 3744, 3848, 3866, 3887, 3888, 3905, 3951, 3954, 4034, 4046, 4066, 4138, 4240, 4256, 4312, 4313, 4334 }; 

 //   _fail_list = {29, 40, 47, 49, 104, 143, 149, 158, 162, 184, 196, 221, 225, 230, 247, 268, 271, 306, 315, 316, 320, 393, 408, 417, 420, 440, 441, 460, 469, 476, 502, 507, 514, 523, 579, 583, 599, 615, 627, 631, 645, 689, 736, 750, 757, 768, 787, 790, 800, 803, 805, 812, 813, 825, 827, 845, 854, 875, 882, 893, 905, 906, 912, 923, 932, 938, 946, 955, 956, 997, 1023, 1026, 1039, 1115, 1116, 1120, 1122, 1138, 1141, 1154, 1193, 1205, 1216, 1223, 1227, 1260, 1299, 1305, 1306, 1313, 1337, 1338, 1355, 1384, 1452, 1461, 1468, 1469, 1477, 1485, 1499, 1529, 1539, 1551, 1554, 1560, 1585, 1588, 1591, 1598, 1623, 1646, 1663, 1665, 1668, 1669, 1672, 1675, 1690, 1691, 1699, 1717, 1735, 1760, 1763, 1764, 1771, 1777, 1786, 1793, 1795, 1798, 1804, 1811, 1820, 1848, 1849, 1866, 1871, 1874, 1878, 1887, 1896, 1903, 1926, 1939, 1946, 1961, 1962, 1966, 2004, 2017, 2028, 2031, 2033, 2058, 2066, 2071, 2074, 2079, 2082, 2083, 2089, 2100, 2167, 2201, 2227, 2239, 2307, 2311, 2314, 2354, 2355, 2375, 2386, 2388, 2393, 2395, 2408, 2427, 2431, 2433, 2440, 2455, 2465, 2478, 2488, 2490, 2492, 2503, 2506, 2513, 2535, 2536, 2550, 2552, 2592, 2631, 2649, 2659, 2663, 2699, 2731, 2736, 2746, 2808, 2830, 2831, 2838, 2851, 2878, 2883, 2889, 2906, 2913, 2914, 2920, 2950, 2955, 2968, 2979, 2980, 2998, 2999, 3011, 3034, 3071, 3086, 3111, 3134, 3138, 3155, 3156, 3165, 3171, 3182, 3184, 3214, 3254, 3277, 3285, 3307, 3338, 3348, 3350, 3352, 3357, 3363, 3373, 3389, 3392, 3398, 3411, 3434, 3444, 3445, 3461, 3465, 3471, 3474, 3484, 3490, 3499, 3514, 3520, 3548, 3558, 3608, 3618, 3620, 3622, 3628, 3640, 3648, 3668, 3681, 3701, 3724, 3731, 3744, 3755, 3774, 3777, 3799, 3848, 3854, 3866, 3869, 3885, 3887, 3888, 3897, 3905, 3925, 3934, 3951, 3954, 3960, 3970, 3988, 4007, 4019, 4032, 4034, 4046, 4051, 4058, 4059, 4063, 4066, 4079, 4109, 4128, 4138, 4147, 4160, 4167, 4177, 4240, 4247, 4256, 4278, 4294, 4304, 4312, 4313, 4316, 4325, 4334, 4338, 4340};

_pi0_list = {29,39,40,47,48,49,59,62,76,79,100,103,104,106,119,143,149,157,158,160,162,171,184,190,196,208,217,221,225,230,247,268,269,271,286,292,306,315,316,319,320,321,364,371,391,393,397,408,417,419,420,424,429,437,440,441,448,452,460,469,476,484,487,489,493,502,504,506,507,514,517,523,533,542,577,579,583,597,599,612,615,620,627,631,645,655,673,681,689,732,736,750,751,757,762,768,776,786,787,790,800,803,805,812,813,825,827,840,845,846,849,853,854,869,875,882,893,897,905,906,912,919,923,931,932,938,943,946,955,956,964,971,997,1023,1026,1034,1037,1039,1045,1047,1090,1093,1115,1116,1120,1122,1123,1129,1138,1141,1144,1147,1154,1161,1190,1193,1205,1216,1219,1223,1225,1227,1243,1260,1275,1285,1299,1305,1306,1308,1313,1325,1328,1337,1338,1355,1359,1384,1385,1391,1394,1404,1411,1430,1438,1446,1452,1461,1468,1469,1475,1477,1485,1490,1499,1509,1516,1529,1533,1539,1541,1545,1551,1554,1560,1561,1573,1585,1588,1591,1598,1612,1623,1643,1645,1646,1651,1661,1663,1665,1666,1668,1669,1671,1672,1675,1690,1691,1699,1707,1709,1717,1726,1735,1736,1760,1763,1764,1771,1777,1780,1786,1793,1795,1798,1801,1804,1805,1811,1820,1845,1848,1849,1859,1861,1866,1871,1874,1878,1886,1887,1896,1897,1903,1919,1926,1939,1940,1944,1945,1946,1949,1961,1962,1966,1983,1984,1995,2004,2017,2027,2028,2031,2033,2037,2044,2049,2056,2058,2066,2071,2074,2079,2082,2083,2089,2100,2121,2122,2132,2141,2151,2163,2167,2171,2174,2195,2201,2215,2223,2227,2236,2239,2259,2263,2265,2285,2307,2311,2314,2320,2338,2352,2354,2355,2372,2375,2386,2388,2391,2393,2395,2398,2408,2413,2417,2420,2424,2427,2428,2431,2433,2435,2440,2455,2465,2471,2478,2488,2489,2490,2492,2496,2498,2503,2506,2513,2535,2536,2550,2552,2560,2573,2576,2592,2604,2615,2631,2649,2659,2663,2698,2699,2709,2731,2736,2739,2746,2799,2808,2830,2831,2838,2851,2875,2878,2883,2889,2897,2906,2913,2914,2920,2933,2950,2955,2968,2979,2980,2981,2986,2995,2998,2999,3000,3011,3019,3031,3034,3057,3071,3075,3077,3078,3086,3111,3118,3123,3127,3134,3135,3138,3143,3155,3156,3164,3165,3171,3182,3184,3204,3206,3207,3214,3225,3246,3249,3254,3265,3267,3277,3285,3306,3307,3312,3327,3338,3346,3348,3350,3352,3357,3362,3363,3373,3381,3389,3392,3398,3411,3414,3433,3434,3444,3445,3461,3465,3471,3474,3484,3490,3499,3504,3514,3517,3520,3542,3547,3548,3552,3558,3581,3601,3608,3611,3612,3618,3620,3622,3628,3638,3640,3648,3668,3681,3701,3718,3724,3731,3740,3744,3755,3774,3777,3783,3784,3799,3827,3828,3844,3846,3848,3854,3860,3861,3866,3869,3880,3885,3887,3888,3897,3898,3905,3908,3918,3925,3927,3934,3949,3951,3954,3959,3960,3963,3970,3974,3984,3988,4007,4019,4032,4034,4046,4051,4058,4059,4062,4063,4066,4079,4090,4104,4109,4128,4135,4136,4138,4141,4147,4148,4157,4160,4167,4177,4184,4205,4226,4229,4240,4247,4256,4263,4278,4294,4304,4312,4313,4316,4325,4334,4338,4340,4341,4352};

    return true;
  }
  
  bool CCpi0Eff::analyze(storage_manager* storage) {

    _event++ ;
    std::cout<<"\nEVENT: "<<_event<<std::endl ;
    //if ( std::find(_pi0_list.begin(),_pi0_list.end(),_event) == _pi0_list.end() )
    //  return false;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
      }

    auto ev_mcv= storage->get_data<event_vertex>("mcvertex"); 
    if(!ev_mcv || !ev_mcv->size() ) {
      std::cout<<"Event has no mcvertex info "<<std::endl;
      return false;
      }

    auto ev_recov= storage->get_data<event_vertex>("numuCC_vertex"); 
    if(!ev_recov || !ev_recov->size() ) {
      std::cout<<"Event has no recovertex info "<<std::endl;
      return false;
      }

    auto rv = ev_recov->at(0); 
    auto mcv = ev_mcv->at(0); 

    auto dist = sqrt( pow(rv.X() - mcv.X(),2) + pow(rv.Y() - mcv.Y(),2) + pow(rv.Z() - mcv.Z(),2));

    auto ev_mcshr= storage->get_data<event_mcshower>("mcreco"); 
    if(!ev_mcshr || !ev_mcshr->size() ) {
      std::cout<<"Event has no mcshower info "<<std::endl;
      return false;
      }

    auto parts = ev_mctruth->at(0).GetParticles();
    bool pi0 = false;
    bool mu  = false ;
    int n_pi0 = 0;
    int n_mes = 0;
    int n_lep = 0;
    int n_mu = 0;
    int n_antimu = 0;
    int n_gamma = 0;

    for ( auto const & p : parts ){
      
      if ( p.StatusCode() == 1 )
        std::cout<<"PDG: "<<p.PdgCode()<<std::endl ;
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){
        pi0 = true;
        n_pi0 ++;
        }   

      if( p.StatusCode() == 1 && p.PdgCode() == 13 ){
        n_mu ++;

	//std::cout<<"Muon Energy:  "<<p.Trajectory().at(0).E()<<", Length: "<<dist<<std::endl;
        mu = true;  
        }   

      if( p.StatusCode() == 1 && p.PdgCode() == -13 )
        n_antimu ++;

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 211 || abs(p.PdgCode()) == 321 ||  p.PdgCode() == 130 
                                   || p.PdgCode() == 310 || abs(p.PdgCode()) == 311 ) ){
        n_mes ++;
	std::cout<<"Meson Energy:  "<<p.Trajectory().at(0).E()<<", Length: "<<dist<<std::endl;
        //break ;
        }   

      if( p.StatusCode() == 1 && (abs(p.PdgCode()) == 11 || p.PdgCode() == -13)){
        n_lep ++;
        //break ;
        }   

      if( p.StatusCode() == 1 && p.PdgCode() == 22)
        n_gamma ++;
      }   

      //std::cout<<"Report: "<<std::endl ;
      //std::cout<<"N pi0: "<<n_pi0<<std::endl;
      //std::cout<<"N meson: "<<n_mes<<std::endl;
      //std::cout<<"N lepton: "<<n_lep<<std::endl;
      //std::cout<<"N mu : "<<n_mu<<std::endl ;

      if( n_mu == 0 && n_pi0 > 0 )
        _nc_pi0 ++; 
      else if( n_mu == 1 && n_pi0 > 0 ){
        _tot_ccpi0 ++; 
      if( n_pi0 > 1)
          _multpi0++;
      else if ( n_mes > 0 ){
          _mesons++;
	  _mes_v.emplace_back(_event);
	  }
      else if ( n_lep > 0)
          _leptons++;
        }
      else if( n_pi0 == 0 && n_gamma > 1 )
        _gammas++;
      //else if ( dist > 10 ) std::cout<<"Event is: "<<_event<<", "<<dist<<std::endl ; 
//      else{
//        std::cout<<"\nEvent is : "<<_event <<", "<<storage->event_id()<<", "<<storage->subrun_id()<<std::endl ;
        //std::cout<<"Report: "<<std::endl ;
        //std::cout<<"N pi0: "<<n_pi0<<std::endl;
        //std::cout<<"N meson: "<<n_mes<<std::endl;
        //std::cout<<"N lepton: "<<n_lep<<std::endl;
        //std::cout<<"N gamma : "<<n_gamma<<std::endl ;
        //std::cout<<"N mu : "<<n_mu<<", "<<n_antimu<<std::endl ;

       //for ( auto const & p : parts )
       //  std::cout<<"PDG: "<<p.PdgCode()<<std::endl ;
//	
//	}

      if( n_mu == 1 && n_pi0 == 1 && n_lep == 0 && n_mes == 0){ 
      //  //std::cout<<"**********************FOUND CCPI0!! "<<std::endl;
      //  _signal++;
        _event_list.emplace_back(_event);
        }   


    return true;
  }

  bool CCpi0Eff::finalize() {

    //std::cout<<"CCpi0 are "<<float(_signal)/(_event)*100<<"\% of BNB ("<<_signal<<"/"<<_event<<")"<<std::endl ;

    std::cout<<"Total CCpi0 : "<<_tot_ccpi0<<"/"<<_fail_list.size()<<std::endl; 

    std::cout<<"NC pi0 : "<<_nc_pi0<<std::endl;
    std::cout<<"Background because mesons: "<<_mesons<<std::endl ;
    std::cout<<"Background because multpi0s: "<<_multpi0<<std::endl ;
    std::cout<<"Background because gammas: "<<_gammas<<std::endl ;
    std::cout<<"Background because leptons: "<<_leptons<<std::endl ;
    std::cout<<"Total accounted backgrounds: "<<_mesons+_multpi0+_gammas+_leptons+_nc_pi0<<std::endl ;

    //std::cout<<_mes_v.size()<<" in Event list :" <<std::endl ;
    //for( auto const & e : _mes_v) std::cout<<e<<", ";

    std::cout<<"\n\n"<<_event_list.size()<<" in Event list :" <<std::endl ;
    for( auto const & e : _event_list) std::cout<<e<<", ";
  
    return true;
  }

}
#endif
