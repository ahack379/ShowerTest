#ifndef LARLITE_GENIEXSECERRORSFULL_CXX
#define LARLITE_GENIEXSECERRORSFULL_CXX

#include "GenieXSecErrorsFull.h"
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

  GenieXSecErrorsFull::GenieXSecErrorsFull() {
   
    _name                    = "GenieXSecErrorsFull";
    _fout                    = 0;
    _tree = nullptr;
    _final_tree = nullptr;

    //_genie_label_v = {"QEMA","QEMA", "NCELaxial", "NCELaxial", "CCResAxial", "CCResAxial", "CCResVector", "CCResVector", "N    CResAxial", "NCResAxial", "NCResVector", "NCResVector", "CohMA", "CohMA", "CohR0", "CohR0", "NonResRvp1pi", "NonResRvp1pi", "NonResRvba    rp1pi", "NonResRvbarp1pi", "NonResRvp2pi", "NonResRvp2pi", "NonResRvbarp2pi", "NonResRvbarp2pi", "ResDecayGamma", "ResDecayGamma", "Res    DecayTheta", "ResDecayTheta", "NC", "NC", "FermiGasModelKf", "FermiGasModelKf", "IntraNukeNmfp", "IntraNukeNmfp", "IntraNukeNcex", "Int    raNukeNcex", "IntraNukeNel", "IntraNukeNel", "IntraNukeNinel", "IntraNukeNinel", "IntraNukeNabs", "IntraNukeNabs", "IntraNukeNpi", "Int    raNukeNpi", "IntraNukePImfp", "IntraNukePImfp", "IntraNukePIcex", "IntraNukePIcex", "IntraNukePIel", "IntraNukePIel", "IntraNukePIinel"    , "IntraNukePIinel", "IntraNukePIabs", "IntraNukePIabs"} ;
    _genie_label_v = {"QEMA", "NCELaxial", "CCResAxial", "CCResVector", "NCResAxial", "NCResVector", "CohMA", "CohR0", "NonResRvp1pi", "NonResRvbarp1pi", "NonResRvp2pi", "NonResRvbarp2pi", "ResDecayGamma", "ResDecayTheta", "NC", "FermiGasModelKf", "IntraNukeNmfp", "IntraNukeNcex", "IntraNukeNel", "IntraNukeNinel", "IntraNukeNabs", "IntraNukeNpi", "IntraNukePImfp", "IntraNukePIcex", "IntraNukePIel", "IntraNukePIinel", "IntraNukePIabs"} ;
    
    fGeometry = nullptr;

  }
  
  bool GenieXSecErrorsFull::initialize(){

    fGeometry = larutil::Geometry::GetME();
    _tot_pot = 0. ;

    _all_evts_nominal = 0;
    _all_evts_m1.resize(30,0) ; // 30 weights
    _all_evts_p1.resize(30,0) ; // 30 weights

    if( !_tree){
       _tree = new TTree("tree","tree");
       _tree->Branch("xsec_mom_truth",&_xsec_mom_truth,"xsec_mom_truth/F");
       _tree->Branch("xsec_theta_truth",&_xsec_theta_truth,"xsec_theta_truth/F");
       _tree->Branch("weight_v",&_weight_v,"std::vector<float>");
       }

    if( !_final_tree){
       _final_tree = new TTree("final_tree","final_tree");
       _final_tree->Branch("all_evts_nominal",&_all_evts_nominal,"all_evts_nominal/F");
       _final_tree->Branch("all_evts_m1",&_all_evts_m1,"std::vector<float>");
       _final_tree->Branch("all_evts_p1",&_all_evts_p1,"std::vector<float>");
    }


    _file.open("list.txt",std::ios_base::in);
 
    return true;
   }

 
  bool GenieXSecErrorsFull::analyze(storage_manager* storage) {


  for(std::string line; std::getline(_file, line); )   //read stream line by line
  {
      std::istringstream in(line);      //make a stream for the line itself
  
      float subrun;
      in >> subrun;                  //and read the first whitespace-separated token
  
      if(subrun == storage->subrun_id() )       //and check its value
      {
          float event ;
          in >> event ;       //now read the whitespace-separated floats
          if ( storage->event_id() == event )
             break;
      }
      else if( storage->subrun_id() > subrun )
        return false; 
   }


    auto ev_pot = storage->get_subrundata<potsummary>("generator"); 
    auto ev_mctruth = storage->get_data<event_mctruth>("generator"); 
    auto ev_mcs = storage->get_data<event_mcshower>("mcreco");
    auto ev_wgt= storage->get_data<event_mceventweight>("genieeventweight"); 

    if( !ev_pot ){ 
       std::cout<<"No POT..."<<std::endl ;
       return false;
       }

    if( storage->subrun_id() != storage->last_subrun_id() )
      _tot_pot += ev_pot->totgoodpot ;

    if(!ev_mctruth || !ev_mctruth->size() ){ 
       std::cout<<"No Truth..."<<std::endl ;
       return false;
     }

    if( !ev_wgt || !ev_wgt->size() ){
      std::cout<<"No event weights..." <<std::endl;
      return false;
     }

    if(!ev_mcs || !ev_mcs->size() ){
      std::cout<<"No mcshower..." <<std::endl;
      return false;
      }
  
    _weight_v.clear();

    auto wgt  = ev_wgt->at(0).GetWeights();

    auto nu  = ev_mctruth->at(0).GetNeutrino();
    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    xyz[0] = traj.at(traj.size() - 1).X();
    xyz[1] = traj.at(traj.size() - 1).Y();
    xyz[2] = traj.at(traj.size() - 1).Z();
    auto nu_energy = traj.at(traj.size() - 1).E();

   if( xyz[0] < 20 || xyz[0] > 236.35 || xyz[1] > 96.5 || xyz[1] < -96.5 || xyz[2] < 10 || xyz[2] > 1026.8 )
        return false;

    auto parts = ev_mctruth->at(0).GetParticles();

    int n_pi0 = 0;
    int n_mu = 0;
    bool isSignal = false ;

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

    // We know in the fv at this point 
    if( n_mu == 1 && n_pi0 == 1 && nu_energy > 0.5 ){
     std::cout<<"WOOOOOOOOOO FOUND ONE "<<std::endl ;

      std::vector<int> shr_ids;
      
      for ( int si = 0; si < ev_mcs->size(); si++){ 

        auto s = ev_mcs->at(si);

        if( s.PdgCode() != 22 ) continue; 
      
        auto st = s.Start();
        auto dist = sqrt( pow(st.X() - start[0],2) + pow(st.Y() - start[1],2) + pow(st.Z() - start[2],2) );
      
        if ( dist < 0.001 )
          shr_ids.emplace_back(si) ;
      }
      
      if( shr_ids.size() == 2 ){

        auto s1 = ev_mcs->at(shr_ids[0]).Start();
        auto s2 = ev_mcs->at(shr_ids[1]).Start();

        auto mag1 = sqrt( s1.Px()*s1.Px()+s1.Py()*s1.Py()+s1.Pz()*s1.Pz() );
        auto mag2 = sqrt( s2.Px()*s2.Px()+s2.Py()*s2.Py()+s2.Pz()*s2.Pz() );
        auto dot = s1.Px()*s2.Px() + s1.Py()*s2.Py() + s1.Pz()*s2.Pz() ;

        lep_dcosz_truth = acos( dot / mag1 / mag2 );  
      }


      _xsec_mom_truth = lep_mom_truth; 
      _xsec_theta_truth = lep_dcosz_truth;
      _all_evts_nominal ++ ;

      auto w_v = wgt.begin()->second; //
      std::cout<<"Number of weights : "<<w_v.size()<<std::endl ;

      for ( int function = 0; function < w_v.size()/2; function++ ){ 

        _all_evts_m1[function] += (w_v.at(2*function)) ; 
        _all_evts_p1[function] += (w_v.at(2*function+1)) ;
        
        _weight_v.emplace_back(w_v.at(2*function));
        _weight_v.emplace_back(w_v.at(2*function+1));
        
        //xsec_mom_truth_p1[function] -> Fill(lep_mom_truth, (m.second.at(2*function+1)));
        //xsec_mom_truth_m1[function] -> Fill(lep_mom_truth, (m.second.at(2*function))); 
        //xsec_theta_truth_p1[function] -> Fill(lep_dcosz_truth, (m.second.at(2*function+1)));
        //xsec_theta_truth_m1[function] -> Fill(lep_dcosz_truth, (m.second.at(2*function)));
       }

      _tree->Fill();

      //xsec_mom_truth -> Fill(lep_mom_truth);
      //xsec_theta_truth -> Fill(lep_dcosz_truth);
      //all_evts_nominal++;
      //for ( int function = 0; function < (m.second.size())/2; function++ ){ 
      //  xsec_mom_truth_p1[function] -> Fill(lep_mom_truth, (m.second.at(2*function+1)));
      //  xsec_mom_truth_m1[function] -> Fill(lep_mom_truth, (m.second.at(2*function))); 
      //  xsec_theta_truth_p1[function] -> Fill(lep_dcosz_truth, (m.second.at(2*function+1)));
      //  xsec_theta_truth_m1[function] -> Fill(lep_dcosz_truth, (m.second.at(2*function)));
      //}

      isSignal = true; 

        }

    return true;
  }

  bool GenieXSecErrorsFull::finalize() {

    std::cout<<"All events: "<<_all_evts_nominal<<std::endl ;
    for( int i = 0 ; i < _all_evts_m1.size(); i++) {
      std::cout<<"\nFunction: "<<_genie_label_v[i]<<std::endl ;
      std::cout<<"All events (-3sig): "<<_all_evts_m1.at(i)<<std::endl ;
      std::cout<<"All events (+3sig): "<<_all_evts_p1.at(i)<<std::endl ;
    }


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
