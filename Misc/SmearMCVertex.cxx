#ifndef LARLITE_SMEARMCVERTEX_CXX
#define LARLITE_SMEARMCVERTEX_CXX

#include "SmearMCVertex.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mctruth.h"

#include "LArUtil/GeometryHelper.h"

#include <random>

namespace larlite {

  bool SmearMCVertex::initialize() {

    _SCE = new larutil::SpaceChargeMicroBooNE();
    _time2cm = larutil::GeometryHelper::GetME()->TimeToCm();;

    return true;
  }
  
  bool SmearMCVertex::analyze(storage_manager* storage) {

    auto ev_mctruth= storage->get_data<event_mctruth>("generator");
    if(!ev_mctruth || !ev_mctruth->size() ) return false;
    
    auto nu = ev_mctruth->at(0).GetNeutrino();
    auto parts = ev_mctruth->at(0).GetParticles();

    auto ev_vtx = storage->get_data<event_vertex>("mcvertex");
    storage->set_id(storage->run_id(), storage->subrun_id(), storage->event_id());
    
    ev_vtx->reserve(1);
    double xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    auto xvtx = traj.at(traj.size() - 1).X();
    auto yvtx = traj.at(traj.size() - 1).Y();
    auto zvtx = traj.at(traj.size() - 1).Z();
    auto tvtx = traj.at(traj.size()-1).T(); // ns
    auto vtxtick = (tvtx / 1000.) * 2.; 
    auto vtxtimecm = vtxtick * _time2cm; 

    // get spacecharge correction
    auto sce_corr = _SCE->GetPosOffsets(xvtx,yvtx,zvtx);
    
    xyz[0] = xvtx + vtxtimecm + _offset - sce_corr.at(0);
    xyz[1] = yvtx + sce_corr.at(1);
    xyz[2] = zvtx + sce_corr.at(2);

    if (_filter) {
      if ( (xyz[0] < 0) || (xyz[0] > 256) || (xyz[1] < -116) || (xyz[1] > 116) || (xyz[2] < 0) || (xyz[2] > 1036) )
        return false;
    }

    // Centering these distributions around 0. std below include 2sig outlier exclusion
    //
    //std::default_random_engine generator;
    std::random_device rd;
    std::mt19937 e2(rd());

    std::normal_distribution<double> x_vtx(0,3.75); //2.77); commented out are the pi0 resolutions
    std::normal_distribution<double> y_vtx(0,3.77); //2.38); comm
    std::normal_distribution<double> z_vtx(0,3.71); //1.97);

    double xnew = x_vtx(e2); //generator);
    double ynew = y_vtx(e2); //generator);
    double znew = z_vtx(e2); //generator);

    //std::cout<<"Vertex: "<<xnew<<", "<<ynew<<", "<<znew<<std::endl ;
    xyz[0] += xnew; 
    xyz[1] += ynew; 
    xyz[2] += znew; 

    
    vertex new_vtx(xyz) ;
    ev_vtx->push_back(new_vtx);

    return true;
  }

  bool SmearMCVertex::finalize() {

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
    // if(_fout) { _fout->cd(); h1->Write(); }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    return true;
  }

}
#endif
