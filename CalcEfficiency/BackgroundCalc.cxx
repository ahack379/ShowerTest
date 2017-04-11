#ifndef LARLITE_CCPI0EFF_CXX
#define LARLITE_CCPI0EFF_CXX

#include "BackgroundCalc.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/shower.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"

namespace larlite {

  bool BackgroundCalc::initialize() {    

    _event = -1; 
    //_event_list.clear();

    _pi0_list = {8,31,35,42,46,48,50,52,55,63,85,90,99,101,104,105,106,113,132,137,146,147,148,155,171,174,181,235,244,250,261,265,276,277,278,294,301,337,370,374,378,397,404,427,430,440,444,449,450,459,473,474,494,500,513,515,537,548,573,579,605,606,614,631,637,650,670,722,734,757,762,798,812,815,840,862,864,898,907,913,953,954,970,971,974,975,980,986,1040,1064,1066,1077,1086,1087,1101,1102,1103,1114,1122,1165,1181,1192,1195,1206,1207,1209,1214,1223,1228,1243,1265,1266,1271,1272,1296,1305,1311,1314,1320,1322,1330,1338,1373,1374,1402,1403,1411,1412,1441,1458,1474,1476,1489,1491,1492,1493,1500,1527,1541,1548,1576,1587,1597,1614,1621,1623,1631,1648,1665,1668,1672,1685,1707,1712,1720,1730,1735,1756,1773,1774,1775,1779,1808,1818,1819,1820,1831,1833,1849,1851,1852,1871,1885,1888,1894,1895,1917,1918,1932,1957,1969,1970,1974,1979,1980,1996,2015,2017,2037,2045,2058,2068,2078,2085,2103,2109,2114,2134,2140,2147,2150,2156,2157,2164,2165,2167,2173,2180,2195,2199,2200,2207,2208,2222,2225,2232,2233,2240,2242,2288,2293,2309,2325,2341,2343,2352,2391,2393,2427,2436,2451,2458,2481,2485,2488,2498,2514,2518,2520,2529,2545,2550,2558,2559,2573,2575,2582,2602,2608,2623,2628,2656,2674,2693,2696,2697,2702,2726,2737,2738,2788,2789,2799,2808,2821,2842,2847,2849,2871,2915,2918,2921,2939,2950,2952,2967,2985,3011,3023,3043,3057,3075,3083,3088,3090,3109,3115,3117,3118,3125,3177,3181,3186,3196,3197,3204,3207,3211,3215,3235,3236,3245,3271,3282,3283,3287,3293,3296,3300,3306,3312,3316,3323,3324,3326,3337,3364,3373,3378,3379,3397,3410,3426,3436,3439,3440,3445,3461,3465,3466,3482,3509,3523,3530,3538,3545,3552,3559,3560,3575,3608,3625,3633,3655,3668,3672,3687,3697,3700,3740,3762,3785,3805,3809,3816,3821,3840,3843,3865,3882,3899,3915,3922,3942,3943,3948,3964,3987,4007,4008,4017,4019,4037,4042,4062,4064,4069,4070,4073,4098,4116,4148,4149,4160,4173,4177,4178,4189,4190,4198,4200,4201,4216,4220,4236,4239,4242,4248,4261,4284,4291,4300,4310,4317,4339,4355,4364,4387,4388,4391,4393,4394,4396,4434,4436,4442,4443,4458,4460,4461,4466,4473,4504,4506,4507,4517,4520,4522,4533,4535,4567,4572,4578,4582,4584,4590,4612,4618,4632,4633,4640,4648,4663,4671,4693,4701,4703,4709,4715,4716,4761,4769,4785,4796,4805,4807,4816,4821,4832,4843,4853,4854,4871,4874,4891,4904,4905,4909,4915,4924,4926,4929,4930,4935,4936,4964,5010,5011,5024,5043,5068,5071,5090,5094,5107,5112,5130,5139,5179,5202,5213,5230,5251,5259,5265,5274,5275,5285,5287,5290,5296,5311,5331,5335,5347,5352,5360,5366,5372,5376,5383,5384,5388,5406,5408,5413,5450,5459,5463,5483,5486,5490,5512,5518,5522,5542,5545,5561,5633,5640,5651,5663,5665,5682,5688,5690,5712,5713,5717,5742,5764,5767,5768,5778,5785,5791,5797,5812,5816,5817,5827,5830,5852,5867,5868,5873,5875,5886,5900,5925,5942,5945,5951,5970,5982,5992,6000,6009,6015,6035,6042,6077,6080,6098,6127,6150,6153,6166,6202,6207,6211,6215,6250,6265,6274,6280,6283,6296,6302,6309,6315,6318,6323,6332,6342,6348,6353,6358,6372,6390,6398,6406,6417,6418,6436,6444,6445,6455,6459,6465,6469,6474,6476,6484,6505,6507,6511,6515,6533,6541,6547,6557,6569,6571,6586,6594,6602,6605,6609,6621,6634,6647,6648,6652,6674,6681,6693,6696,6762,6768,6775,6776,6780,6809,6818,6829,6830,6833,6836,6843,6848,6852,6883,6897,6899,6926,6927,6953,6954,6955,6957,6964,6967,6973,6974,6983,6993,7004,7035,7037,7047,7058,7060,7062,7074,7107,7111,7132,7148,7162,7177,7184,7187,7195,7209,7214,7215,7233,7239,7248,7258,7262,7286,7293,7308,7309,7311,7323,7332,7361,7363,7391,7395,7402,7407,7408,7411,7412,7424,7431,7441,7447,7452,7460,7463,7473,7479,7481,7512,7517,7521,7527,7535,7536,7546,7555,7564,7583,7603,7609,7619,7641,7645,7669,7680,7683,7693,7715,7725,7731,7732,7739,7751,7755,7756,7764,7772,7775,7792,7794,7810,7812,7826,7828,7829,7833,7836,7837,7843,7867,7874,7883,7894,7906,7917,7924,7936,7942,7978,8001,8006,8024,8025,8066,8072,8079,8080,8091,8102};

    _n_noise = 0;     // 1
    _n_cosmic = 0;    // 2
    _n_nue = 0;       // 3
    _n_antinumu = 0;  // 4
    _n_nc = 0;        // 5
    _multpi0 = 0;     // 6
    _ccpi0_outfv = 0; // 7 
    _tot_ccpi0 = 0;   // 8
    _n_gammas = 0;    // 9
    _n_ccother = 0;   // 10

    if ( !_tree ){
      _tree = new TTree("tree","");
      _tree->Branch("event",&_event,"event/I");
      _tree->Branch("bkgd_id",&_bkgd_id,"bkgd_id/I");

      _tree->Branch("vtx_x",&_vtx_x,"vtx_x/F");
      _tree->Branch("vtx_y",&_vtx_y,"vtx_y/F");
      _tree->Branch("vtx_z",&_vtx_z,"vtx_z/F");
      _tree->Branch("mu_angle",&_mu_angle,"mu_angle/F");
      _tree->Branch("mu_len",&_mu_len,"mu_len/F");

      _tree->Branch("pi0_mass",&_pi0_mass,"pi0_mass/F");
      _tree->Branch("pi0_oangle",&_pi0_oangle,"pi0_oangle/F");
      _tree->Branch("pi0_mom",&_pi0_mom,"pi0_mom/F");
      _tree->Branch("pi0_low_shrE",&_pi0_low_shrE,"pi0_low_shrE/F");
      _tree->Branch("pi0_high_shrE",&_pi0_high_shrE,"pi0_high_shrE/F");
      _tree->Branch("pi0_low_radL",&_pi0_low_radL,"pi0_low_radL/F");
      _tree->Branch("pi0_high_radL",&_pi0_high_radL,"pi0_high_radL/F");
   }


  //std::cout<<"PI0 LIST!  "<<_pi0_list.size()<<std::endl;

    return true;
  }
  
  bool BackgroundCalc::analyze(storage_manager* storage) {

    _bkgd_id = -1 ;
    _event++ ;
    if ( std::find(_pi0_list.begin(),_pi0_list.end(),_event) == _pi0_list.end() )
      return false;

    //std::cout<<"\nEvent is : "<<_event <<", "<<storage->event_id()<<", "<<storage->subrun_id()<<std::endl ;

    auto ev_mctruth= storage->get_data<event_mctruth>("generator"); 
    if(!ev_mctruth || !ev_mctruth->size() ) {
      std::cout<<"Event has no mctruth info "<<std::endl;
      return false;
      }

    auto ev_mctrk = storage->get_data<event_mctrack>("mcreco");
    if ( !ev_mctrk || !ev_mctrk->size() ) {std::cout<<"No MCTrack!" <<std::endl ; return false; }

    auto ev_vtx= storage->get_data<event_vertex>("numuCC_vertex"); 
    if(!ev_vtx || !ev_vtx->size() ) {
      std::cout<<"Event has no recovertex info "<<std::endl;
      return false;
      }

    auto ev_trk = storage->get_data<event_track>("pandoraNu");
    if ( !ev_trk || !ev_trk->size() ) {std::cout<<"No Track!" <<std::endl ; return false; }

    auto ev_tagged_trk = storage->get_data<event_track>("numuCC_track");
    if ( !ev_tagged_trk || !ev_tagged_trk->size() ){ std::cout<<"No Tagged Track!" <<std::endl ; return false; }

    auto ev_s = storage->get_data<event_shower>("showerreco");

    if( !ev_s || !ev_s->size() || ev_s->size() < 2 ){
      std::cout<<"Not enough reco'd showers..." <<std::endl;
      return false;
     }   

    for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

        auto const& shr1 = ev_s->at(s1);

        for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

            if (s2 <= s1) continue;

            auto const& shr2 = ev_s->at(s2);

            geoalgo::Vector_t rev_shr1(-1.*shr1.Direction()) ;
            geoalgo::Vector_t rev_shr2(-1.*shr2.Direction()) ;

            // Make the backwards projection for the showers
            auto shr1_bkwrd_hl = ::geoalgo::HalfLine_t(shr1.ShowerStart(),rev_shr1);
            auto shr2_bkwrd_hl = ::geoalgo::HalfLine_t(shr2.ShowerStart(),rev_shr2);

            // Calc the Opening angle of the showers
            double oangle = acos( shr1.Direction().Dot(shr2.Direction())) ;

            // Calc the vertex point of the two showers. the true designated backwards project
            geoalgo::Point_t vertex(3);

            auto st1 = shr1.ShowerStart();
            auto st2= shr2.ShowerStart();
            auto dir1 = shr1.Direction();
            auto dir2 = shr2.Direction();
            geoalgo::HalfLine_t shr1_hl(st1.X(),st1.Y(),st1.Z(),dir1.X(),dir1.Y(), dir1.Z() );
            geoalgo::HalfLine_t shr2_hl(st2.X(),st2.Y(),st2.Z(),dir2.X(),dir2.Y(), dir2.Z() );

            _geoAlgo.commonOrigin(shr1_hl, shr2_hl, vertex, true);

            // Calc Diretion of two correlated shower
            geoalgo::Vector_t momentum(3);// need to fill out
            geoalgo::Vector_t mom_vect(shr2.Direction()*shr1.Energy(2) +shr1.Direction()*shr2.Energy(2)) ;
            auto tot_pi0_mom = sqrt(pow(mom_vect[0],2) + pow(mom_vect[1],2) + pow(mom_vect[2],2) );

            //===========================================
            auto radL_shr1 = vertex.Dist(shr1.ShowerStart());
            auto radL_shr2 = vertex.Dist(shr2.ShowerStart());
            //===========================================

            if( oangle < 0.35 ) continue;
            if( pow( _geoAlgo.SqDist(shr1_bkwrd_hl, shr2_bkwrd_hl), 0.5 ) > 4.) continue;
            if( radL_shr1 > 62. || radL_shr2 > 62. ) continue; 

            _pi0_mass      = sqrt(2 * shr1.Energy() * shr2.Energy() *(1.-cos(oangle)));
            _pi0_mom       = tot_pi0_mom;
            _pi0_oangle    = oangle;
            _pi0_low_shrE  = shr1.Energy() < shr2.Energy() ? shr1.Energy() : shr2.Energy() ;
            _pi0_high_shrE = shr1.Energy() < shr2.Energy() ? shr2.Energy() : shr1.Energy() ;
            _pi0_low_radL  = shr1.Energy() < shr2.Energy() ? radL_shr1 : radL_shr2 ;
            _pi0_high_radL = shr1.Energy() < shr2.Energy() ? radL_shr2 : radL_shr1 ;
        }// shower ID 2 
      }// shower ID 1 

      _mu_angle      = cos(ev_tagged_trk->at(0).Theta());
      _mu_len        = ev_tagged_trk->at(0).Length(0); // Calculates the length from point 0 to end

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Want to be able to access the origin of the tagged muon. Thus, need to find it, and 
    // Ask for its origin.  Need to match to MCtrack to do this
    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();

    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();

    auto vtx = ev_vtx->at(0); 
    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };
    _vtx_x = vtx.X();
    _vtx_y = vtx.Y();
    _vtx_z = vtx.Z();

    auto vtx_diff = sqrt(pow(xyz[0] - _vtx_x,2) + pow(xyz[1] - _vtx_y,2) + pow(xyz[2] - _vtx_z,2));

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

       //if( _event == 1103 && (dist_st < 6 || dist_end < 6 ) )
       //   std::cout<<"1103 : "<<dist_st<<", "<<dist_end <<std::endl ;

       if ( dist_st < 3 || dist_end < 3 ){

          float len = ev_trk->at(ti).Length();
          trk_map.emplace(1./len,ti);
	  min_dist = dist_st < dist_end ? dist_st : dist_end ; 
        }   
     }   

    int max_it = -1;
    float max_dot = -1.;

    if( trk_map.size() ) { 

    auto m_st = ev_tagged_trk->at(0).VertexDirection();     
    auto m_norm = sqrt( pow(m_st.Px(),2) + pow(m_st.Py(),2) + pow(m_st.Pz(),2)); 

      for( auto & ti : trk_map ){
            
        auto t = ev_trk->at(ti.second);
        auto t_st = t.VertexDirection();
        
        auto dot = (m_st.Px() * t_st.Px() + m_st.Py() * t_st.Py() + m_st.Pz() * t_st.Pz())/m_norm ;

        if ( fabs(dot) > max_dot ){
             max_dot = dot;
             max_it = ti.second ;
          }
        }
      }   
    else return false ;

    // Annnnnnd same thing for MC, to grab the origin of the track and assess backgrounds properly
    std::multimap<float,int> mctrk_map ;
    auto tag_trk = ev_tagged_trk->at(0);
    auto tag_st = tag_trk.Vertex() ;
    auto tag_end = tag_trk.End() ;
    float mc_min_dist = 1e9;

    for ( size_t ti = 0; ti < ev_mctrk->size(); ti++ ) { 

      auto mc_vtx = ev_mctrk->at(ti).Start() ;
      auto mc_end = ev_mctrk->at(ti).End() ;
    
      float dist_st = sqrt(  pow(mc_vtx.X() - tag_st.X(),2) + 
                             pow(mc_vtx.Y() - tag_st.Y(),2) + 
                             pow(mc_vtx.Z() - tag_st.Z(),2) );  

      float dist_end = sqrt( pow(mc_vtx.X() - tag_end.X(),2) + 
                             pow(mc_vtx.Y() - tag_end.Y(),2) + 
                             pow(mc_vtx.Z() - tag_end.Z(),2) );  

       //if ( (dist_st < 100 || dist_end < 100) && _event == 2208  ) 
       //  std::cout<<"dist_st end: "<<dist_st<<", "<<dist_end<<", "<<tag_end.Y()<<", "<<tag_st.Y()<<std::endl ;

       if ( dist_st < 16 || dist_end < 16 ){
          float length = sqrt( pow(mc_end.X() - mc_vtx.X(),2) + 
                           pow(mc_end.Y() - mc_vtx.Y(),2) + 
                           pow(mc_end.Z() - mc_vtx.Z(),2) );  
          mctrk_map.emplace(1./length,ti);
	  mc_min_dist = dist_st < dist_end ? dist_st : dist_end ; 
        }   
     }   
     float temp = mctrk_map.begin()->first; 


    int mc_max_it = -1;
    float mc_max_dot = -1.;

    if( mctrk_map.size() ) { 

    auto tag_st = tag_trk.VertexDirection();     
    auto tag_norm = sqrt( pow(tag_st.Px(),2) + pow(tag_st.Py(),2) + pow(tag_st.Pz(),2)); 

      for( auto & ti : mctrk_map ){
            
        auto mc = ev_mctrk->at(ti.second);
        auto mc_st = mc.Start();
        auto mc_norm = sqrt( pow(mc_st.Px(),2) + pow(mc_st.Py(),2) + pow(mc_st.Pz(),2) );
        
        auto dot = (tag_st.Px() * mc_st.Px() + tag_st.Py() * mc_st.Py() + tag_st.Pz() * mc_st.Pz())/tag_norm / mc_norm ;

        if ( fabs(dot) > mc_max_dot ){
             mc_max_dot = dot;
             mc_max_it = ti.second ;
          }
        }
      }   
    // If no true tracks aligned with reco track, mark it as noise
    else {
       std::cout<<"\nEvent is : "<<_event <<", mult: "<<trk_map.size()<<std::endl ; //storage->event_id()<<", "<<storage->subrun_id()<<std::endl ;
       //std::cout<<"Vertex difference: "<<vtx_diff<<std::endl ;
       //std::cout<<"MC length : "<<1./temp <<", and ntracks: "<<ev_mctrk->size()<<std::endl ;

       //if( trk_map.size() >= 2 ){

       //   auto it0 = trk_map.begin()->second ;
       //   auto it1_0 = ++trk_map.begin();

       //   auto it1 = it1_0->second ;
       //   
       //   auto t0 = ev_trk->at(it0);
       //   auto t1 = ev_trk->at(it1);

       //   auto d0 = t0.VertexDirection();
       //   auto d1 = t1.VertexDirection();

       //   auto dot = d0.Px() * d1.Px() + d0.Py() * d1.Py() + d0.Pz() * d1.Pz() ;
       //   std::cout<<"DOT: "<<dot <<", "<<t0.Vertex().X()<<", "<<t0.Vertex().Y()<<", "<<t0.Vertex().Z()<<", "
       //            <<t0.End().X()<<", "<<t0.End().Y()<<", "<<t0.End().Z()<<std::endl ; 
       //   std::cout<<"LENGTH: "<<t0.Length()<<", "<<t1.Length()<<std::endl ; 
       //}

      //_n_noise++;
      _n_cosmic++;
      _bkgd_id = 2 ;
      _tree->Fill();

     return false;

     //auto parts = ev_mctruth->at(0).GetParticles();
     // for ( auto const & p : parts )
     //   if ( p.StatusCode() == 1 && p.PdgCode() < 1000 ) std::cout<<"PDG: "<<p.PdgCode()<<std::endl ;
    }

      auto mc_vtx = ev_mctrk->at(mc_max_it).Start() ;
      auto mc_end = ev_mctrk->at(mc_max_it).End() ;
      float length = sqrt( pow(mc_end.X() - mc_vtx.X(),2) + 
                           pow(mc_end.Y() - mc_vtx.Y(),2) + 
                           pow(mc_end.Z() - mc_vtx.Z(),2) );  

      //std::cout<<"MC CHOSEN TRK LENGHT : "<<length<<std::endl ;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Now have only events with true tagged tracks.  Next check containment of vertex
    ///////////////////////////////////////////////////////////////////////////////////////////////


    bool infv = true;
    //auto dist = sqrt( pow(vtx.X() - xyz[0],2) + pow(vtx.Y() - xyz[1],2) + pow(vtx.Z() - xyz[2],2));

    if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 106.5 || xyz[1] < -106.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
      infv = false;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /// Now count the number of backgrounds and signals
    ///////////////////////////////////////////////////////////////////////////////////////////////
    auto parts = ev_mctruth->at(0).GetParticles();
    int n_pi0 = 0;
    int n_mu = 0;
    int n_gamma = 0;

    // It turns out that for MC files, origin 2 is cosmic, 1 is neutrino. Can check this using 
    // BNB Only files or OpenCosmic files. The opposite is the case for reco tracks through back tracker, for some
    // reason.
    if( ev_mctrk->at(mc_max_it).Origin() == 2 ){
      //std::cout<<"\nEvent is : "<<_event <<", "<<mc_max_dot<<std::endl ; //storage->event_id()<<", "<<storage->subrun_id()<<std::endl ;
      _n_cosmic++;
      _bkgd_id = 2; 
    }
    else if( abs(nu.Nu().PdgCode()) == 12){
      _n_nue ++ ;
      _bkgd_id = 3;
    } 
    else if( nu.Nu().PdgCode() == -14){
      _n_antinumu++ ;
      _bkgd_id = 4;
    }
    else if( nu.Nu().PdgCode() == 14 && nu.CCNC() == 1 ){
      _n_nc++;
      _bkgd_id = 5;
    }

    if( _bkgd_id == -1 ){
      for ( auto const & p : parts ){
     
        if( p.StatusCode() == 1 && p.PdgCode() == 111 )
          n_pi0 ++;
        if( p.StatusCode() == 1 && p.PdgCode() == 13 )
          n_mu ++;
        if( p.StatusCode() == 1 && p.PdgCode() == 22)
          n_gamma ++;

        }   

        if( n_mu == 1){
          if( n_pi0 > 1 ){
            _bkgd_id = 6; 
            _multpi0++;
          }
          else if( n_pi0 == 1 && !infv ){
            _bkgd_id = 7;
            _ccpi0_outfv ++ ;
          }
          else if( n_pi0 == 1 && infv ){ 
            _bkgd_id = 8;
            _tot_ccpi0 ++; 
            _event_list.emplace_back(_event);
            }
          else if( n_pi0 == 0 && n_gamma > 1 ){
            _bkgd_id = 9;
            _n_gammas++;
          }
          else{
            _bkgd_id = 10;
            _n_ccother ++;
	    _ccother_list.emplace_back(_event);
          }
        }
      }
    
    _tree->Fill();    

    return true;
  }

  bool BackgroundCalc::finalize() {

    std::cout<<"Signals: "<<std::endl ;
    std::cout<<"Total CCpi0 : "<<_tot_ccpi0<<"/"<<_event_list.size()<<std::endl; 

    // Note that cc other includes secondary pi0s.
    std::cout<<"\nBackgrounds: "<<std::endl;
    std::cout<<"1) Noise : "<<_n_noise<< std::endl;
    std::cout<<"2) BNB Cosmic : "<<_n_cosmic<< std::endl;
    std::cout<<"3) Nue : "<<_n_nue<<std::endl;
    std::cout<<"4) Antinumus: "<<_n_antinumu<<std::endl;
    std::cout<<"5) NC pi0 : "<<_n_nc<<std::endl; // nc_pi0<<std::endl;
    std::cout<<"6) Multpi0s: "<<_multpi0<<std::endl ;
    std::cout<<"7) Nu outsude FV: "<<_ccpi0_outfv <<std::endl ;
    std::cout<<"8) Ngamma Events: "<<_n_gammas<<std::endl ;
    std::cout<<"9) OTHER CC Events: "<<_n_ccother<<std::endl ;

    std::cout<<"Total accounted backgrounds: "<< _n_noise + _n_cosmic + _n_nue + _n_antinumu + _n_nc +
             _multpi0 + _ccpi0_outfv + _n_gammas + _n_ccother <<std::endl ;

    //std::cout<<"\n\n"<<_event_list.size()<<" in Event list :" <<std::endl ;
    //for( auto const & e : _event_list) std::cout<<e<<", ";

    std::cout<<"\n\n"<<_ccother_list.size()<<" in CCOther list :" <<std::endl ;
    for( auto const & e : _ccother_list) std::cout<<e<<", ";

    if ( _fout ){
      _fout->cd();
      _tree->Write();
    }
  
    return true;
  }

}
#endif
