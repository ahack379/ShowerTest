#ifndef LARLITE_SHOWERCALIB_CXX
#define LARLITE_SHOWERCALIB_CXX

#include "ShowerCalib.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/shower.h"
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/simch.h"
#include "DataFormat/event_ass.h"
#include "LArUtil/DetectorProperties.h"
#include "LArUtil/Geometry.h"

#include <algorithm>

namespace larlite {

  bool ShowerCalib::initialize() {
    
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
      _gamma_tree->Branch("true_adj_gamma_e",&_true_adj_gamma_e,"true_adj_gamma_pi0_e/F");
      _gamma_tree->Branch("reco_adj_gamma_e",&_reco_adj_gamma_e,"reco_adj_gamma_pi0_e/F");
      _gamma_tree->Branch("true_reco_dot",&_true_reco_dot,"true_reco_dot/F");
      _gamma_tree->Branch("true_rad_l",&_true_rad_l,"true_rad_l/F");
      _gamma_tree->Branch("reco_rad_l",&_reco_rad_l,"reco_rad_l/F");
      }

    _event = -1;

  _event_list = {69 }; 
  
   std::cout<<"Event list size: "<<_event_list.size() <<std::endl ;

    return true;
  }

  void ShowerCalib::Clear(){
    _true_pi0_e = -999;
    _true_angle = -9;
    _true_asym = -9;
    _reco_pi0_e = -999;

    _n_true_pi0 = 0;

    _true_gamma_e = -9;
    _reco_gamma_e = -9;
    _true_adj_gamma_e = -9;
    _reco_adj_gamma_e = -9;
    _true_rad_l = -9;
    _reco_rad_l = -9;
    _true_reco_dot = -9;
    }
  
  bool ShowerCalib::analyze(storage_manager* storage) {

    _event ++; 
    Clear() ;

    //if(std::find(_event_list.begin(), _event_list.end(), _event) == _event_list.end()) 
    //  return false;

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

  // Get clusters associated to shower
  auto ev_hit = storage->get_data<event_hit>("hit02");
  if ( !ev_hit || !ev_hit->size() ) return false;

  auto const& ev_clus = storage->get_data<event_cluster>("ImageClusterHit");

  auto ev_ass_c = storage->get_data<larlite::event_ass>("ImageClusterHit");
  auto const& ass_hit_v = ev_ass_c->association(ev_clus->id(), ev_hit->id());

  auto ev_s = storage->get_data<event_shower>("showerreco"); 
  if( !ev_s || !ev_s->size() ){
    std::cout<<"No shower..." <<std::endl;
    return false;
    }

  if( ev_s->size() < 2 ) return false ;

  auto ev_ass = storage->get_data<larlite::event_ass>("showerreco");
  auto const& ass_clus_v = ev_ass->association(ev_s->id(), ev_clus->id());

  auto geom    = ::larutil::Geometry::GetME();
  auto ev_simch   = storage->get_data<event_simch>("largeant");

  std::map<unsigned int, size_t> Ch2SimchIdx_map;
  for (size_t i=0; i < ev_simch->size(); i++)
    Ch2SimchIdx_map[ ev_simch->at(i).Channel() ] = i;

  std::map<int,float> shrid_to_energyloss_reco ;
  std::map<int,float> shrid_to_energyloss_true;

  // Loop over showers. I think.
  for (size_t i = 0; i < ass_clus_v.size(); i++ ){
      
    for (size_t j = 0; j < ass_clus_v.at(i).size(); j++ ){
        auto clus_id = ass_clus_v.at(i).at(j); 
        auto iclus = ev_clus->at(clus_id);
       
        if (iclus.Plane().Plane != 2 ) continue; 

        auto const & ass_clus_hit = ass_hit_v.at(clus_id) ;

        std::cout<<"Number of hits: "<<ass_clus_hit.size()<<std::endl;

        float sub_true_charge = 0 ;
        float sub_reco_charge = 0 ;

        for (auto const& k : ass_clus_hit){

          auto const& reco_hit = ev_hit->at ( k );

          float tick = reco_hit.PeakTime();
          float wire = reco_hit.WireID().Wire;
          float reco_area = reco_hit.Integral();

          auto clocktick = larutil::DetectorProperties::GetME()->SamplingRate() * 1.e-3; //time sample in usec
          float lifetime_corr = exp( tick * clocktick / 8000);
          float electrons = reco_area * 189 ; //mcc8 value

          float dQ = electrons * lifetime_corr * 23.6 * 10e-6 ; 
          float dE = dQ / 0.62 ; // 0.62 -> recomb factor

          // find the simchannel info for this channel
          auto chan = geom->PlaneWireToChannel(iclus.Plane().Plane,wire);

          // are there simch @ this channel?
          if (Ch2SimchIdx_map.find(chan) == Ch2SimchIdx_map.end())
            continue;

          // find the simch associated with this vector of LArLite IDEs
          auto const& simch = ev_simch->at( Ch2SimchIdx_map[chan] );

          // get the charge deposited in the appropriate time-range
          auto const& ide_v = simch.TrackIDsAndEnergies( reco_hit.PeakTime() - 4*reco_hit.RMS() + 7298,
                                                         reco_hit.PeakTime() + 4*reco_hit.RMS() + 7298);

          auto q = 0;
          float temp_energy = 0;
          for (auto const&  ide : ide_v){
            q += ide.numElectrons;
            temp_energy +=  ide.energy ;
              }

          auto gain = float(q) / reco_area ; 

          if ( gain > 240 ){
            std::cout<<"GAIN IS HUGE: "<<gain<<", "<<sub_true_charge<<", "<<sub_reco_charge<<std::endl ;
	    sub_true_charge += temp_energy;
            sub_reco_charge += dE;
             } 
           
             }
          shrid_to_energyloss_true[i] = sub_true_charge ;
          shrid_to_energyloss_reco[i] = sub_reco_charge ;
           }
        }

    //std::cout<<"True and reco charge: "<<sub_true_charge<<", "<<sub_reco_charge<<std::endl ;
 
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
           _reco_adj_gamma_e = _reco_gamma_e - shrid_to_energyloss_reco[0] ; 
           _true_adj_gamma_e = _true_gamma_e - shrid_to_energyloss_true[0]; 
           std::cout<<"Subtracting off energy loss: "<<shrid_to_energyloss_reco[0]<<std::endl  ;

	   _true_reco_dot = max_dot;
	   _gamma_tree->Fill();
	   std::cout<<"Reco and tru gamma e: "<<_reco_gamma_e<<", "<<_true_gamma_e<<std::endl ;
	   std::cout<<"Start of shower: "<<mcs_1.Start().X()<<", "<<mcs_1.Start().Y() <<", "<<mcs_1.Start().Z()<<std::endl ;
	   std::cout<<"End of shower: "<<mcs_1.End().X()<<", "<<mcs_1.End().Y() <<", "<<mcs_1.End().Z()<<std::endl ;

	   auto mcs_2 = ev_mcs->at(min_mcs) ;
	   auto recos_2 = ev_s->at(min_recos) ;
	   auto dot = mcs_2.Start().Px() * recos_2.Direction().Px() + mcs_2.Start().Py() * recos_2.Direction().Py() + mcs_2.Start().Pz() * recos_2.Direction().Pz() ;
	   auto mag_reco = sqrt( pow(recos_2.Direction().Px(),2) + pow(recos_2.Direction().Py(),2) + pow(recos_2.Direction().Pz(),2) ); 
	   auto mag_mcs = sqrt( pow(mcs_2.Start().Px(),2) + pow(mcs_2.Start().Py(),2) + pow(mcs_2.Start().Pz(),2) ); 
	   dot /= (mag_mcs * mag_reco);

	   _true_gamma_e = mcs_2.Start().E();
	   _reco_gamma_e = recos_2.Energy(2);

           _reco_adj_gamma_e = _reco_gamma_e - shrid_to_energyloss_reco[1] ; 
           _true_adj_gamma_e = _true_gamma_e - shrid_to_energyloss_true[1]; 

	   _true_reco_dot = dot;
	   _gamma_tree->Fill();

	   std::cout<<"\nReco and tru gamma e: "<<_reco_gamma_e<<", "<<_true_gamma_e<<std::endl ;
           std::cout<<"Subtracting off energy loss: "<<shrid_to_energyloss_reco[1]<<std::endl  ;
	   std::cout<<"Start of shower: "<<mcs_2.Start().X()<<", "<<mcs_2.Start().Y() <<", "<<mcs_2.Start().Z()<<std::endl ;
	   std::cout<<"End of shower: "<<mcs_2.End().X()<<", "<<mcs_2.End().Y() <<", "<<mcs_2.End().Z()<<std::endl ;
         }
  
    return true;
  }

  bool ShowerCalib::finalize() {

    if(_fout) { 
      _fout->cd(); 
      _pi0_tree->Write(); 
      _gamma_tree->Write(); 
      }
  
    return true;
  }

}
#endif
