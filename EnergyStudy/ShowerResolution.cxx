#ifndef LARLITE_SHOWERRESOLUTION_CXX
#define LARLITE_SHOWERRESOLUTION_CXX

#include "ShowerResolution.h"
#include "LArUtil/GeometryHelper.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  bool ShowerResolution::initialize() {

   if ( !_tree ){
      _tree = new TTree("tree","tree");
      _tree->Branch("e_mcc",&_e_mcc,"e_mcc/F");
      _tree->Branch("e_true",&_e_true,"e_true/F");
    }   

   _event = -1 ;

    return true;
  }
  
  bool ShowerResolution::analyze(storage_manager* storage) {

    _event++; 
    std::cout<<"\nEVENT: "<<_event <<std::endl ;

    auto const& geomH = ::larutil::GeometryHelper::GetME();

    //auto ev_mcc = storage->get_data<event_cluster>("mccluster") ;
    //if( !ev_mcc || !ev_mcc->size() ) return false; 

    auto ev_mcs = storage->get_data<event_mcshower>("mcreco") ;
    if( !ev_mcs || !ev_mcs->size() ) return false; 

    auto ev_ass = storage->get_data<larlite::event_ass>("mccluster");
    auto const& ass_keys = ev_ass->association_keys();

    if ( ass_keys.size() == 0 ){  return false;}

    larlite::event_hit *ev_hit = nullptr;
    auto ass_mcclus_v = storage->find_one_ass( ass_keys[0].second, ev_hit, ev_ass->name() );

    if (!ev_hit || ev_hit->size() == 0){ 
      std::cout << "No ass! exit" << std::endl;
      return false;
    }   

    if (ass_mcclus_v.size() == 0){ 
      std::cout << "No ass! exit" << std::endl;
      return false;
    }   

    auto ev_mctruth = storage->get_data<event_mctruth>("generator") ;
    if(!ev_mctruth || !ev_mctruth->size() ) { 
      std::cout<<"Event has no mctruth info "<<std::endl; 
      return false;
    }   

    auto & truth = ev_mctruth->at(0);
    auto & nu  = truth.GetNeutrino();
    auto traj = nu.Nu().Trajectory();
    auto xvtx = traj.at(traj.size() - 1).X();
    auto yvtx = traj.at(traj.size() - 1).Y();
    auto zvtx = traj.at(traj.size() - 1).Z();

    // Now deal with shower comparisons
    std::vector<int> shr_ids;

    for ( int si = 0; si < ev_mcs->size(); si++){

      auto s = ev_mcs->at(si);
      auto st = s.Start();
      auto dist = sqrt( pow(st.X() - xvtx,2) + pow(st.Y() - yvtx,2) + pow(st.Z() - zvtx,2) );

      if ( dist < 0.0001 && s.DetProfile().E() > 0 && s.MotherPdgCode() == 111 ){ shr_ids.emplace_back(si) ; } 
    }   

    if ( shr_ids.size() != 2 ) {std::cout<<" shr size: "<<shr_ids.size() <<std::endl ; return false; }

    std::vector<std::pair<int,float>> clus_e_v ;
    std::map<int,float> clus_m ;

    int temp_it = 0;

    // Fill map with hits from mccluster : clusterID
    for( int i = 0; i < ass_mcclus_v.size(); i ++ ){

     bool wrongplane = false ;
     float clus_e = 0. ;

     for ( int j = 0; j < ass_mcclus_v[i].size(); j++ ){

       auto h = ev_hit->at(ass_mcclus_v[i][j]);

       if ( h.WireID().Plane != 2){ 
         wrongplane = true;
         break ;
       }

       auto clocktick = larutil::DetectorProperties::GetME()->SamplingRate() * 1.e-3; //time sample in usec

       float lifetime_corr = exp( h.PeakTime() * clocktick / 1.e20);
       float electrons = h.Integral() * 198.; //mcc8 value
       float dQ = electrons * lifetime_corr * 23.6 * 1e-6 ;
       float dE = dQ / 0.577 ; // 0.62 -> recomb factor

       clus_e += dE ;

     } 

     if ( wrongplane ) continue ;

     clus_e_v.emplace_back(std::make_pair(temp_it,clus_e)) ;
     temp_it ++ ;
     }

    if ( clus_e_v.size() != 2 ) return false ;
     
    if ( clus_e_v.at(0).second > clus_e_v.at(1).second ){
      clus_m[0] = clus_e_v.at(0).second ; 
      clus_m[1] = clus_e_v.at(1).second ; 
    }
    else {
      clus_m[0] = clus_e_v.at(1).second ; 
      clus_m[1] = clus_e_v.at(0).second ; 
    }

    

    auto s0 = ev_mcs->at(shr_ids[0]).DetProfile().E();
    auto s1 = ev_mcs->at(shr_ids[1]).DetProfile().E() ;

     std::map<int,float> e_m ;

    if ( s0 > s1 ){
      e_m[0] = s0 ;
      e_m[1] = s1 ;
    }
    else {
      e_m[0] = s1 ;
      e_m[1] = s0 ;
    }

      _e_mcc = clus_m[1] ; //e_m[1] ;
      _e_true = e_m[1] ;

      std::cout<<"Energy : "<<_e_mcc<<", "<<_e_true<<std::endl; 
      _tree->Fill();

      _e_mcc = clus_m[0] ;
      _e_true = e_m[0] ;

      std::cout<<"Energy : "<<_e_mcc<<", "<<_e_true<<std::endl; 

      _tree->Fill();



    //int max_s_it = -1; 
    //int max_c_it = -1; 
    //int max_dot = -10; 

    //for ( int i = 0; i < shr_ids->size(); i++ ){

    //  auto s = ev_mcs->at(shr_ids[si]);

    //  if( s.PdgCode() != 22 ) continue; 
    //
    //  auto s_i = ev_mcs->at(0); 
    //  std::vector<double> st = { s_i.Start().X(), s_i.Start().Y(), s_i.Start().Z() };

    //  auto stWT  = geomH->Point_3Dto2D(st,2);
    //  auto st_w = stWT.w ; // Comes in cm; 
    //  auto st_t = stWT.t + 800. * geomH->TimeToCm() ; // 800 for single particle files

    //  std::vector<double> end = { s_i.End().X(), s_i.End().Y(), s_i.End().Z() };

    //  auto endWT  = geomH->Point_3Dto2D(end,2);
    //  auto end_w = endWT.w ; // Comes in cm; 
    //  auto end_t = endWT.t + 800. * geomH->TimeToCm() ; // 800 for single particle files

    //  auto dist_mc = sqrt( pow( end_w - st_w,2 ) + pow( end_t - st_t,2) );
    //  std::vector<float> mc_dir= { (end_w - st_w) / dist_mc , (end_t - st_t) / dist_mc } ;
    //  
    //  for ( int j = 0; j < ev_mcc->size(); j++ ){
    //    auto c_j = ev_mcc->at(j);

    //    if ( c_j.PlaneID().Plane != 2 ) continue;

    //    auto cend_w = c_j.EndWire();
    //    auto cend_t = c_j.EndTick();

    //    auto cst_w = c_j.StartWire();
    //    auto cst_t = c_j.StartTick();

    //    auto dist = sqrt( pow( cend_w - cst_w,2 ) + pow( cend_t - cst_t,2) );

    //    std::vector<float> dir = { (cend_w - cst_w) / dist , (cend_t - cst_t) / dist } ;

    //    auto dot = dir.at(0) * dir_mc.at(0) + dir.at(1) * dir_mc.at(1) ;

    //    if ( dot > max_dot ){
    //      max_dot = dot ;
    //      max_s_it = shr_ids[si]; 
    //      max_c_it = c_j ;

    //     }

    //  }
    //}

    //auto min_c_it = max_c_it == 0 ? 1 : 0 ;
    //auto min_s_it = max_s_it == 0 ? 1 : 0 ;

 

  
    return true;
  }

  bool ShowerResolution::finalize() {

    if(_fout) { _fout->cd(); _tree->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    return true;
  }

}
#endif
