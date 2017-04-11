#ifndef LARLITE_NUEVARIABLESTUDY_CXX
#define LARLITE_NUEVARIABLESTUDY_CXX

#include "NueVariableStudy.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/vertex.h"
#include "DataFormat/potsummary.h"
#include <algorithm>

namespace larlite {

  bool NueVariableStudy::initialize() {
    
    if( !_nue_selection ){
      _nue_selection = new TTree("_nue_selection","");
      _nue_selection->Branch("_energy",&_energy,"_energy/F");
      _nue_selection->Branch("_mc_xvtx",&_mc_xvtx,"_mc_xvtx/F");
      _nue_selection->Branch("_mc_yvtx",&_mc_yvtx,"_mc_yvtx/F");
      _nue_selection->Branch("_mc_zvtx",&_mc_zvtx,"_mc_zvtx/F");
      _nue_selection->Branch("_rc_xvtx",&_rc_xvtx,"_rc_xvtx/F");
      _nue_selection->Branch("_rc_yvtx",&_rc_yvtx,"_rc_yvtx/F");
      _nue_selection->Branch("_rc_zvtx",&_rc_zvtx,"_rc_zvtx/F");
      _nue_selection->Branch("_vtx_dist",&_vtx_dist,"_vtx_dist/F");
      _nue_selection->Branch("_pass",&_pass,"_pass/I");
      _nue_selection->Branch("_nc",&_nc,"_nc/B");
    }

    if( !_pot_tree){
      _pot_tree = new TTree("_pot_tree","");
      _pot_tree->Branch("_pottotl",&_pottotl,"_pottotl/F");
      _pot_tree->Branch("_evttotl",&_evttotl,"_evttotl/F");
      }

    _event = -1;
    _pottotl = 0;

    _event_list = {3,66,82,89,111,117,133,135,161,166,169,175,176,180,199,203,211,217,218,220,225,228,232,242,251,262,277,289,300,304,305,312,320,321,324,327,346,357,385,391,400,410,411,416,420,432,440,447,448,457,478,493,508,534,547,551,587,589,606,618,639,644,648,657,660,681,701,708,728,746,749,763,769,771,790,791,794,821,824,839,841,842,860,868,880,889,940,948,952,965,975,997,1001,1048,1058,1065,1088,1099,1108,1111,1113,1116,1120,1124,1126,1134,1140,1144,1147,1159,1177,1182,1184,1206,1210,1223,1224,1245,1247,1248,1265,1306,1311,1317,1318,1321,1354,1365,1366,1367,1371,1382,1385,1390,1406,1407,1409,1411,1427,1435,1444,1446,1452,1454,1477,1495,1520,1545,1561,1571,1606,1613,1614,1622,1635,1647,1663,1665,1688,1690,1716,1718,1732,1735,1738,1758,1773,1786,1804,1813,1815,1817,1824,1825,1840,1841,1850,1852,1863,1889,1901,1911,1912,1925,1935,1945,1956,1959,1960,1965,1976,1984,1993,1997,2011,2020,2036,2050,2053,2060,2063,2068,2070,2073,2085,2086,2091,2103,2112,2124,2143,2158,2174,2180,2191,2202,2204,2210,2213,2216,2217,2224,2226,2229,2236,2249,2257,2266,2274,2283,2285,2295,2299,2309,2311,2317,2318,2327,2343,2355,2359,2371,2380,2387,2393,2394,2397,2403,2410,2457,2459,2484,2511,2515,2518,2529,2540,2542,2554,2558,2564,2582,2593,2605,2610,2614,2617,2642,2651,2662,2665,2682,2689,2696,2701,2708,2712,2726,2729,2748,2755,2759,2772,2791,2813,2820,2824,2838,2842,2850,2879,2882,2900,2915,2923,2925,2930,2932,2944,2957,2964,2967,2977,2990,3003,3015,3018,3024,3029};

    return true;
  }

  bool NueVariableStudy::analyze(storage_manager* storage) {

    _event++; 
    std::cout<<"\nEvent is : "<<_event <<std::endl ;

    auto ev_pot = storage->get_subrundata<potsummary>("generator"); 
    if( storage->subrun_id() != storage->last_subrun_id() )
      _pottotl += ev_pot->totgoodpot ;

    if(std::find(_event_list.begin(), _event_list.end(), _event) == _event_list.end()) 
      return false;


    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ){
      std::cout<<"No mctruth..." <<std::endl;
      return false;
      }

    auto ev_recovtx = storage->get_data<event_vertex>("numuCC_vertex"); 
    if(!ev_recovtx || !ev_recovtx->size() ){
      std::cout<<"No recovtx..." <<std::endl;
      return false;
      }

    auto & nu  = ev_mctruth->at(0).GetNeutrino();

    if( abs(nu.Nu().PdgCode()) != 12 ) return false;

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    if( xyz[0] < 0 || xyz[0] > 256.35 || xyz[1] < -116.5 || xyz[1] > 116.5 || xyz[2] < 0 || xyz[2] > 1036.8 )
      return false ;
    
    _mc_xvtx = xyz[0];
    _mc_yvtx = xyz[1];
    _mc_zvtx = xyz[2];
    _energy = ev_mctruth->at(0).GetNeutrino().Nu().Trajectory().at(0).E() ;
    _nc = nu.CCNC() == 1 ? 0 : 1;

    auto rc_vtx = ev_recovtx->at(0);
    _rc_xvtx = rc_vtx.X();
    _rc_yvtx = rc_vtx.Y();
    _rc_zvtx = rc_vtx.Z();

    _vtx_dist = sqrt( pow(_mc_xvtx - _rc_xvtx,2) + pow(_mc_yvtx - _rc_yvtx,2) + pow(_mc_zvtx - _rc_zvtx,2));

    _pass = 10 ;
     
    _nue_selection->Fill();
      
    return true;
  }

  bool NueVariableStudy::finalize() {

    _evttotl = _event+1 ;
    _pot_tree->Fill();

    if(_fout) { 
      _fout->cd(); 
      _nue_selection->Write();
      _pot_tree->Write();
      }
  
    return true;
  }

}
#endif
