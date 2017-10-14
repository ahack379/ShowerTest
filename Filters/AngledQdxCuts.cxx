#ifndef LARLITE_ANGLEDQDXCUTS_CXX
#define LARLITE_ANGLEDQDXCUTS_CXX

#include "AngledQdxCuts.h"
#include "TMath.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/calorimetry.h"

namespace larlite {

  bool AngledQdxCuts::initialize() {

    return true;
  }

  double AngledQdxCuts::Median(std::vector<double> input){

    int N = input.size();

    double median;

    std::sort(input.begin(), input.end());
    if (N % 2 == 0){ median = (input[((N/2) - 1)] + input[N/2]) / 2;}
    else if( N == 1){median = input[N];}
    else{            median = input[N/2];}
    return median;
  }

  double AngledQdxCuts::TrunMean(std::vector <double> poop){

    double RMS = TMath::RMS(poop.begin(),poop.end());
    double median = Median(poop);

    std::vector<double> TLMean;

    for(int i = 0; i < int(poop.size()); i++){
      if(poop[i] < median+RMS && poop[i] > median-RMS){TLMean.push_back(poop[i]);}
    }

    return TMath::Mean(TLMean.begin(), TLMean.end());
  }

  double AngledQdxCuts::MaxDeflection(::larlite::track trk){

    double max = 0;
    std::vector<TVector3> mom;

    for(int i = 0; i < trk.NumberTrajectoryPoints(); i++){

        TVector3 nextpos;
        TVector3 nextmom;
        trk.TrajectoryAtPoint(i,nextpos,nextmom);
        mom.push_back(nextmom);
      
    }

    for(int i = 0; i < int(mom.size())-1; i++){
      if(max < mom.at(i).Angle(mom.at(i+1))) max = mom.at(i).Angle(mom.at(i+1));
    }

    return max*(180./3.14159265);

  }
  
  bool AngledQdxCuts::analyze(storage_manager* storage) {


      auto ev_t_p = storage->get_data<event_track>("pandoraNu");

      auto ev_t = storage->get_data<event_track>("numuCC_track");
      auto trk = ev_t->at(0);
      auto tag_st = trk.Vertex() ;

      float min_dist = 1e9;
      int min_it = -1;

      for ( int i = 0; i < ev_t_p->size(); i++){

        auto t = ev_t_p->at(i);
        auto st = t.Vertex() ;
        auto dist = sqrt( pow(st.X() - tag_st.X(),2) + pow(st.Y() - tag_st.Y(),2) + pow(st.Z() - tag_st.Z(),2) );
	if ( dist < min_dist ){
	  min_dist = dist;
	  min_it = i;
	}
      }

    auto ev_calo= storage->get_data<event_calorimetry>("pandoraNucalo"); 

    if ( !ev_calo || ev_calo->size() == 0 ) {
      std::cout << "No such calo associated to track! " << std::endl;
      return false;
    }   

    auto ev_ass = storage->get_data<larlite::event_ass>("pandoraNucalo"); 

    if ( !ev_ass || ev_ass->size() == 0 ) {
      std::cout << "No such association! " << std::endl;
      return false;
    }

    auto const& ass_calo_v = ev_ass->association(ev_t_p->id(), ev_calo->id());
    if ( ass_calo_v.size() == 0) {
      std::cout << "No ass from track => hit! " << std::endl;
      return false;
    }   

    auto TrackMaxDeflection = MaxDeflection(trk);

    if ( TrackMaxDeflection > 8 ){
      std::cout<<"Failed track deflection: "<<TrackMaxDeflection<<std::endl;
      return false;
    }
    
    auto len = trk.Length();
    
    int N = 0;
    std::vector<double> dqdx; 

    // Get calo ID for plane 2 at the tagged track.
    auto calo_it = ass_calo_v.at(min_it).at(2) ;
    auto calo_i = ev_calo->at(calo_it); 

    for(int i = 0; i < calo_i.dQdx().size(); i++){

      if ( calo_i.dQdx().at(i) <= 0 ) continue;
      N++;

      if ( _use_mcbnbcos)
        dqdx.push_back(calo_i.dQdx().at(i) * 198.);
      else 
        dqdx.push_back(calo_i.dQdx().at(i) * 243.);
    }
    
      //for( int i = 0; i < int(trk.NumberdQdx(geo::kW)); i++){   
      //
      //  if(trk.DQdxAtPoint(i,geo::kW) <= 0) continue; 
      //  N++;

      //  if ( _use_mcbnbcos)
      //    dqdx.push_back(trk.DQdxAtPoint(i,geo::kW) * 198.);
      //  else 
      //    dqdx.push_back(trk.DQdxAtPoint(i,geo::kW) * 243.);
      //
      //}   

      if(N == 0){ 
        dqdx.clear(); 
        std::cout<<"No Points!"<<std::endl;
        ev_t->clear();
	return false; 
      }

      std::sort(dqdx.begin(),dqdx.end());

      auto TrackTLMeandQdx = TrunMean(dqdx);
      dqdx.clear();

      if( len > 40 && TrackTLMeandQdx < 70000. )
        return true;
      else{
        std::cout<<"FAILED DQDX "<<N<<", "<<len<<", "<<TrackTLMeandQdx<<std::endl;
        ev_t->clear();
        return false;
      }

    }

  bool AngledQdxCuts::finalize() {

    return true;
  }

}
#endif
