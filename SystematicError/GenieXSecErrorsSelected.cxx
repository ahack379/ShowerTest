#ifndef LARLITE_GENIEXSECERRORSSELECTED_CXX
#define LARLITE_GENIEXSECERRORSSELECTED_CXX

#include "GenieXSecErrorsSelected.h"
#include "DataFormat/track.h"
#include "DataFormat/opflash.h"
#include "DataFormat/vertex.h"
#include "DataFormat/calorimetry.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/hit.h"
#include "DataFormat/shower.h"
#include "DataFormat/potsummary.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mceventweight.h"

// HitRatio
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

namespace larlite {

  GenieXSecErrorsSelected::GenieXSecErrorsSelected() {
   
    _name                    = "GenieXSecErrorsSelected";
    _fout                    = 0;
    _tree = nullptr;
    _final_tree = nullptr;

//    _genie_label_v = {"QEMA", "NCELaxial", "CCResAxial", "CCResVector", "NCResAxial", "NCResVector", "CohMA", "CohR0", "NonResRvp1pi", "NonResRvbarp1pi", "NonResRvp2pi", "NonResRvbarp2pi", "ResDecayGamma", "ResDecayTheta", "NC", "FermiGasModelKf", "IntraNukeNmfp", "IntraNukeNcex", "IntraNukeNel", "IntraNukeNinel", "IntraNukeNabs", "IntraNukeNpi", "IntraNukePImfp", "IntraNukePIcex", "IntraNukePIel", "IntraNukePIinel", "IntraNukePIabs"} ;

   _genie_label_v = {"FermiGasModelKf", "FermiGasModelSf", "IntraNukeNabs", "IntraNukeNcex", "IntraNukeNel", "IntraNukeNinel", "IntraNukeNmfp", "IntraNukeNpi", "IntraNukePIabs", "IntraNukePIcex", "IntraNukePIel", "IntraNukePIinel", "IntraNukePImfp", "IntraNukePIpi", "NC", "NonResRvbarp1pi", "NonResRvbarppi", "NonResRvp1pi", "NonResRvppi", "ResDecayEta", "ResDecayGamma", "ResDecayTheta", "ccresAxial", "ccresVector", "cohMA", "cohR0", "ncelAxial", "ncelEta", "ncresAxial", "ncresVector", "qema", "qevec"}; 
    
    fGeometry = nullptr;

    _events = 0;
  }
  
  bool GenieXSecErrorsSelected::initialize(){

    fGeometry = larutil::Geometry::GetME();
    _tot_pot = 0. ;
    int funcs = 32; //64 total, +- for each func

    _sel_evts_nominal = 0;
    _sel_evts_m1.resize(funcs,0) ; // 30 weights
    _sel_evts_p1.resize(funcs,0) ; // 30 weights

    _bkgd_evts_nominal = 0;
    _bkgd_evts_m1.resize(funcs,0) ; // 30 weights
    _bkgd_evts_p1.resize(funcs,0) ; // 30 weights

    if( !_tree){
       _tree = new TTree("tree","tree");
       _tree->Branch("xsec_mom_truth",&_xsec_mom_truth,"xsec_mom_truth/F");
       _tree->Branch("xsec_theta_truth",&_xsec_theta_truth,"xsec_theta_truth/F");
       _tree->Branch("weight_v",&_weight_v,"std::vector<float>");
       }

    if( !_final_tree){
       _final_tree = new TTree("final_tree","final_tree");
       _final_tree->Branch("sel_evts_nominal",&_sel_evts_nominal,"sel_evts_nominal/F");
       _final_tree->Branch("sel_evts_m1",&_sel_evts_m1,"std::vector<float>");
       _final_tree->Branch("sel_evts_p1",&_sel_evts_p1,"std::vector<float>");
       _final_tree->Branch("bkgd_evts_nominal",&_bkgd_evts_nominal,"bkgd_evts_nominal/F");
       _final_tree->Branch("bkgd_evts_m1",&_bkgd_evts_m1,"std::vector<float>");
       _final_tree->Branch("bkgd_evts_p1",&_bkgd_evts_p1,"std::vector<float>");
    }

    return true;
   }

 
  bool GenieXSecErrorsSelected::analyze(storage_manager* storage) {

    auto ev_mctruth = storage->get_data<event_mctruth>("generator"); 
    //auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    auto ev_wgt= storage->get_data<event_mceventweight>("genieeventweight"); 

    if(!ev_mctruth || !ev_mctruth->size() ){ 
       std::cout<<"No Truth..."<<std::endl ;
       return false;
     }

    if( !ev_wgt || !ev_wgt->size() ){
      std::cout<<"No event weights..." <<std::endl;
      return false;
     }

    _events ++ ;


    //if(!ev_mcs || !ev_mcs->size() ){
    //  std::cout<<"No mcshower..." <<std::endl;
    //  return false;
    //  }
  
    auto nu  = ev_mctruth->at(0).GetNeutrino();
    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();
    auto nu_energy = traj.at(traj.size() - 1).E();

    bool infv = true;

   if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
      infv = false; 

    auto parts = ev_mctruth->at(0).GetParticles();

    int n_pi0 = 0;
    int n_mu = 0;

    float lep_mom_truth = 0., lep_dcosz_truth = 0. ;
    std::vector<float> start(3,0) ;
    
    for ( auto const & p : parts ){
    
      if( p.StatusCode() == 1 && p.PdgCode() == 111 ){
        n_pi0 += 1;
        lep_mom_truth = sqrt( pow(p.Trajectory().at(0).Px(),2) + pow(p.Trajectory().at(0).Py(),2) + 
	                      pow(p.Trajectory().at(0).Pz(),2) )*1000;
        start[0] = p.Trajectory().at(0).X();
        start[1] = p.Trajectory().at(0).Y();
        start[2] = p.Trajectory().at(0).Z();
        }

      if( p.StatusCode() == 1 && p.PdgCode() == 13 )
        n_mu += 1;
    }

    _weight_v.clear();
    auto wgt  = ev_wgt->at(0).GetWeights();
//    _events++ ;

    // We know in the fv at this point 
    if( n_mu == 1 && n_pi0 == 1 && nu_energy > 0.5 && infv){

      //std::vector<int> shr_ids;
      //for ( int si = 0; si < ev_mcs->size(); si++){ 

      //  auto s = ev_mcs->at(si);
      //  if( s.PdgCode() != 22 ) continue; 
      //
      //  auto st = s.Start();
      //  auto dist = sqrt( pow(st.X() - start[0],2) + pow(st.Y() - start[1],2) + pow(st.Z() - start[2],2) );
      //
      //  if ( dist < 0.001 )
      //    shr_ids.emplace_back(si) ;
      //}
      //
      //if( shr_ids.size() == 2 ){

      //  auto s1 = ev_mcs->at(shr_ids[0]).Start();
      //  auto s2 = ev_mcs->at(shr_ids[1]).Start();
      //  auto mag1 = sqrt( s1.Px()*s1.Px()+s1.Py()*s1.Py()+s1.Pz()*s1.Pz() );
      //  auto mag2 = sqrt( s2.Px()*s2.Px()+s2.Py()*s2.Py()+s2.Pz()*s2.Pz() );
      //  auto dot = s1.Px()*s2.Px() + s1.Py()*s2.Py() + s1.Pz()*s2.Pz() ;
      //  lep_dcosz_truth = acos( dot / mag1 / mag2 );  
      //}

      _xsec_mom_truth = lep_mom_truth; 
      _xsec_theta_truth = lep_dcosz_truth;
      _sel_evts_nominal ++ ;

      //auto w_v = wgt.begin()->second; //
      //std::cout<<"Number of weights : "<<w_v.size()<<std::endl ;

      //for ( int function = 0; function < w_v.size()/2; function++ ){ 

      //  _sel_evts_m1[function] += (w_v.at(2*function)) ; 
      //  _sel_evts_p1[function] += (w_v.at(2*function+1)) ;
      //  _weight_v.emplace_back(w_v.at(2*function));
      //  _weight_v.emplace_back(w_v.at(2*function+1));
      //  
      // }

      int it = 0;
      for ( auto const & m : wgt ) { 
         auto w_v = m.second ;
         std::cout<<"Parameter: "<<m.first<<", "<<m.second.size() <<std::endl;
         //for ( auto const & w : m.second){
           _sel_evts_m1[it] += (w_v.at(0)) ; 
           _sel_evts_p1[it] += (w_v.at(1)) ;
         it++;
         // }   
       }

      _tree->Fill();

      }
      else{
        _bkgd_evts_nominal ++ ;

        int it = 0;
        for ( auto const & m : wgt ) { 
           auto w_v = m.second ;
             _bkgd_evts_m1[it] += (w_v.at(0)) ; 
             _bkgd_evts_p1[it] += (w_v.at(1)) ;
           it++;
         }   

        //auto w_v = wgt.begin()->second; //
        //for ( int function = 0; function < w_v.size()/2; function++ ){ 
        //  _bkgd_evts_m1[function] += (w_v.at(2*function)) ; 
        //  _bkgd_evts_p1[function] += (w_v.at(2*function+1)) ;
        //  
        // }
      }

    return true;
  }

  bool GenieXSecErrorsSelected::finalize() {

    std::cout<<"All events: "<<_sel_evts_nominal<<", "<<_bkgd_evts_nominal<<std::endl ;

    std::cout<<"Events: "<<_events <<std::endl ;
    //for( int i = 0 ; i < _sel_evts_m1.size(); i++) {
    //  std::cout<<"\nFunction: "<<_genie_label_v[i]<<std::endl ;
    //  std::cout<<"All events (-3sig): "<<_sel_evts_m1.at(i)<<std::endl ;
    //  std::cout<<"All events (+3sig): "<<_sel_evts_p1.at(i)<<std::endl ;
    //}

    for( int i = 0 ; i < _sel_evts_m1.size(); i++) 
      std::cout<<_sel_evts_m1[i]<<", " ;

    std::cout<<std::endl ;
    for( int i = 0 ; i < _sel_evts_m1.size(); i++) 
      std::cout<<_sel_evts_p1[i]<<", " ;

    std::cout<<std::endl ;

    std::cout<<"BACKGROUND: "<<std::endl ;
    for( int i = 0 ; i < _sel_evts_m1.size(); i++) 
      std::cout<<_bkgd_evts_m1[i]<<", " ;

    std::cout<<std::endl ;
    for( int i = 0 ; i < _sel_evts_m1.size(); i++) 
      std::cout<<_bkgd_evts_p1[i]<<", " ;

    std::cout<<std::endl ;

    _final_tree->Fill();
    if(_fout){
     _fout->cd();
     _tree->Write();
     _final_tree->Write();
    }
     
    return true;
  }

}
#endif
