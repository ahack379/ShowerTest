#ifndef LARLITE_HITTHRESHSTUDY_CXX
#define LARLITE_HITTHRESHSTUDY_CXX

#include "HitThreshStudy.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/track.h"
#include "DataFormat/wire.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/simch.h"
#include "DataFormat/rawdigit.h"
#include "DataFormat/vertex.h"
#include "LArUtil/Geometry.h"

namespace larlite {

  bool HitThreshStudy::initialize() {

    _event = -1;

    _tree = new TTree("tree","Signal Processing Hit Normalization");
    _tree->Branch("_pl",&_pl,"pl/I");
    _tree->Branch("_wire",&_wire,"wire/I");
    _tree->Branch("_chan",&_chan,"chan/I");
    _tree->Branch("_reco_area",&_reco_area,"reco_area/D");
    _tree->Branch("_reco_ampl",&_reco_ampl,"reco_ampl/D");
    _tree->Branch("_q",&_q,"q/D");
    _tree->Branch("_tick",&_tick,"tick/I");
    _tree->Branch("_event",&_event,"event/I");
    _tree->Branch("_clus",&_clus,"clus/I");
    _tree->Branch("_x_st",&_x_st,"x_st/D");
    _tree->Branch("_y_st",&_y_st,"y_st/D");
    _tree->Branch("_z_st",&_z_st,"z_st/D");
    _tree->Branch("_px_st",&_px_st,"px_st/D");
    _tree->Branch("_py_st",&_py_st,"py_st/D");
    _tree->Branch("_pz_st",&_pz_st,"pz_st/D");
    _tree->Branch("_e_st",&_e_st,"e_st/D");
    _tree->Branch("_x_dp",&_x_dp,"x_dp/D");
    _tree->Branch("_y_dp",&_y_dp,"y_dp/D");
    _tree->Branch("_z_dp",&_z_dp,"z_dp/D");
    _tree->Branch("_px_dp",&_px_dp,"px_dp/D");
    _tree->Branch("_py_dp",&_py_dp,"py_dp/D");
    _tree->Branch("_pz_dp",&_pz_dp,"pz_dp/D");
    _tree->Branch("_e_dp",&_e_dp,"e_dp/D");

    return true;
  }
  
  bool HitThreshStudy::analyze(storage_manager* storage) {

    _event ++ ;
    std::cout<<"Event: "<<_event<<std::endl ;

    auto geom    = ::larutil::Geometry::GetME();

    auto const& ev_mcs = storage->get_data<event_mcshower>("mcreco");

    if ( !ev_mcs or (ev_mcs->size() == 0) ){
      std::cout << "no mcs found. exiting" << std::endl;
      return false;
    }

    _x_st = ev_mcs->at(0).Start().X();
    _y_st = ev_mcs->at(0).Start().Y();
    _z_st = ev_mcs->at(0).Start().Z();
    _px_st = ev_mcs->at(0).Start().Px();
    _py_st = ev_mcs->at(0).Start().Py();
    _pz_st = ev_mcs->at(0).Start().Pz();
    _e_st = ev_mcs->at(0).Start().E();

    _x_dp = ev_mcs->at(0).DetProfile().X();
    _y_dp = ev_mcs->at(0).DetProfile().Y();
    _z_dp = ev_mcs->at(0).DetProfile().Z();
    _px_dp = ev_mcs->at(0).DetProfile().Px();
    _py_dp = ev_mcs->at(0).DetProfile().Py();
    _pz_dp = ev_mcs->at(0).DetProfile().Pz();
    _e_dp = ev_mcs->at(0).DetProfile().E();
   
    //auto const& ev_clus = storage->get_data<event_cluster>(_clus_producer);
    //if ( !ev_clus or (ev_clus->size() == 0) ){
    //  std::cout << "no hits found. exiting" << std::endl;
    //  return false;
    //}
    //event_hit* ev_recohit = nullptr;
    //auto const& ass_hit_v = storage->find_one_ass(ev_clus->id(), ev_recohit, ev_clus->name());
    //
    //if ( !ev_recohit or (ev_recohit->size() == 0) ){
    //  std::cout << "no hits found. exiting" << std::endl;
    //  return false;
    //}
    
    auto const& ev_hit = storage->get_data<event_hit>("gaushit");

    if ( !ev_hit or (ev_hit->size() == 0) ){
      std::cout << "no hits found. exiting" << std::endl;
      return false;
    }
    

    // load simchannels
    auto ev_simch   = storage->get_data<event_simch>("largeant");
    // make a map from channel to position in ev_simch
    std::map<unsigned int, size_t> Ch2SimchIdx_map;
    for (size_t i=0; i < ev_simch->size(); i++)
      Ch2SimchIdx_map[ ev_simch->at(i).Channel() ] = i;
    
    // make a map from reco hit index to wire number
    std::map<size_t, int> _HitIdx_to_WireNum;
    // make a map from wire number to vector of raw hit indices
    std::map<int, std::vector<size_t> > _WireNum_to_RawHitIdx_v;

    _reco_area = 0;
    _q = 0;
   
    //for (size_t i = 0; i < ass_hit_v.size(); i++ ){

    //  _clus = i ;
    //  auto const & ass_clus_hit = ass_hit_v.at(i) ;

    //  if (ass_clus_hit.size() < 20)
    //    continue;
    //  
    //  //std::cout<<"Number of hits: "<<ass_clus_hit.size()<<std::endl;

    //for (auto const& i : ass_clus_hit){
    for (size_t i = 0; i < ev_hit->size(); i++ ){
	
	//auto const& reco_hit = ev_recohit->at ( i );
	auto const& reco_hit = ev_hit->at ( i );
	if ( reco_hit.WireID().Plane != 2 ) continue ;
	
	_tick = reco_hit.PeakTime();
	_wire = reco_hit.WireID().Wire;
	_pl   = reco_hit.WireID().Plane;
	_reco_ampl = reco_hit.PeakAmplitude();
	_reco_area += reco_hit.Integral();

	if ( _reco_area <= 0.01 ) continue;
	
	// find the simchannel info for this channel
	auto chan = geom->PlaneWireToChannel(_pl,_wire);

	_chan = chan ;
	
	// are there simch @ this channel?
	if (Ch2SimchIdx_map.find(chan) == Ch2SimchIdx_map.end())
	  continue;
	
	// find the simch associated with this 
	// vector of LArLite IDEs
	auto const& simch = ev_simch->at( Ch2SimchIdx_map[chan] );
	
	// get the charge deposited in the appropriate time-range
	//auto const& ide_v = simch.TrackIDsAndEnergies( 7298 + 700, 7298 + 7000 );
	auto const& ide_v = simch.TrackIDsAndEnergies( reco_hit.PeakTime() - 2*reco_hit.RMS() + 7298, reco_hit.PeakTime() + 2*reco_hit.RMS() + 7298);
	
	for (auto const&  ide : ide_v)
	  _q += ide.numElectrons;

      }// for all reconstructed hits
    //}// for all clusters

    _tree->Fill();

    return true;
  }

  bool HitThreshStudy::finalize() {

    _tree->Write();

    return true;
  }

}
#endif
