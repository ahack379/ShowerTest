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

  std::vector<float> x_res={8.0 , 2.0 , 6.0 , 13.0 , 6.0 , 8.0 , 10.0 , 11.0 , 8.0 , 4.0 , 9.0 , 13.0 , 9.0 , 13.0 , 12.0 , 6.0 , 17.0 , 13.0 , 10.0 , 8.0 , 9.0 , 10.0 , 9.0 , 5.0 , 14.0 , 9.0 , 10.0 , 12.0 , 11.0 , 8.0 , 11.0 , 14.0 , 9.0 , 11.0 , 13.0 , 15.0 , 13.0 , 9.0 , 12.0 , 8.0 , 15.0 , 8.0 , 18.0 , 19.0 , 17.0 , 11.0 , 22.0 , 14.0 , 18.0 , 17.0 , 15.0 , 15.0 , 17.0 , 24.0 , 14.0 , 14.0 , 19.0 , 16.0 , 17.0 , 27.0 , 29.0 , 32.0 , 30.0 , 32.0 , 37.0 , 41.0 , 60.0 , 80.0 , 98.0 , 234.0 , 725.0 , 2413.0 , 6176.0 , 14288.0 , 5258.0 , 1859.0 , 536.0 , 179.0 , 93.0 , 69.0 , 45.0 , 29.0 , 34.0 , 36.0 , 34.0 , 21.0 , 27.0 , 16.0 , 19.0 , 19.0 , 22.0 , 20.0 , 15.0 , 12.0 , 18.0 , 16.0 , 22.0 , 15.0 , 13.0 , 9.0 , 13.0 , 14.0 , 14.0 , 7.0 , 9.0 , 7.0 , 19.0 , 17.0 , 8.0 , 10.0 , 11.0 , 12.0 , 8.0 , 12.0 , 14.0 , 14.0 , 16.0 , 16.0 , 7.0 , 8.0 , 12.0 , 10.0 , 9.0 , 17.0 , 13.0 , 8.0 , 10.0 , 6.0 , 11.0 , 14.0 , 7.0 , 3.0 , 8.0 , 9.0 , 7.0 , 16.0 , 7.0 , 8.0 , 11.0 , 8.0 , 14.0 , 13.0 , 10.0 , 5.0 , 8.0 , 8.0 , 14.0 , 9.0 , 12.0 , 12.0};

  std::vector<float> y_res = {8.0 , 3.0 , 8.0 , 5.0 , 10.0 , 11.0 , 7.0 , 5.0 , 13.0 , 5.0 , 9.0 , 6.0 , 7.0 , 5.0 , 5.0 , 13.0 , 8.0 , 12.0 , 15.0 , 13.0 , 13.0 , 9.0 , 14.0 , 14.0 , 11.0 , 12.0 , 10.0 , 18.0 , 8.0 , 17.0 , 13.0 , 17.0 , 18.0 , 16.0 , 17.0 , 12.0 , 11.0 , 12.0 , 22.0 , 18.0 , 15.0 , 12.0 , 15.0 , 20.0 , 20.0 , 14.0 , 18.0 , 21.0 , 26.0 , 19.0 , 21.0 , 18.0 , 25.0 , 21.0 , 25.0 , 21.0 , 24.0 , 46.0 , 23.0 , 30.0 , 28.0 , 28.0 , 33.0 , 41.0 , 29.0 , 58.0 , 65.0 , 67.0 , 92.0 , 147.0 , 285.0 , 555.0 , 1226.0 , 3280.0 , 10983.0 , 10932.0 , 2787.0 , 1021.0 , 416.0 , 234.0 , 134.0 , 85.0 , 66.0 , 51.0 , 37.0 , 40.0 , 52.0 , 36.0 , 28.0 , 25.0 , 36.0 , 22.0 , 31.0 , 23.0 , 19.0 , 23.0 , 18.0 , 22.0 , 17.0 , 20.0 , 14.0 , 19.0 , 15.0 , 12.0 , 16.0 , 17.0 , 23.0 , 10.0 , 13.0 , 14.0 , 8.0 , 18.0 , 13.0 , 22.0 , 12.0 , 13.0 , 9.0 , 13.0 , 12.0 , 18.0 , 9.0 , 9.0 , 11.0 , 13.0 , 12.0 , 14.0 , 8.0 , 11.0 , 7.0 , 9.0 , 5.0 , 10.0 , 7.0 , 13.0 , 9.0 , 16.0 , 8.0 , 5.0 , 10.0 , 11.0 , 10.0 , 7.0 , 6.0 , 6.0 , 11.0 , 3.0 , 8.0 , 6.0 , 11.0 , 7.0};

  std::vector<float> z_res = {5.0 , 9.0 , 3.0 , 7.0 , 5.0 , 9.0 , 9.0 , 2.0 , 4.0 , 5.0 , 0.0 , 6.0 , 3.0 , 7.0 , 9.0 , 6.0 , 7.0 , 6.0 , 6.0 , 4.0 , 9.0 , 6.0 , 5.0 , 2.0 , 6.0 , 5.0 , 4.0 , 9.0 , 5.0 , 11.0 , 8.0 , 3.0 , 6.0 , 4.0 , 5.0 , 5.0 , 5.0 , 5.0 , 7.0 , 9.0 , 7.0 , 3.0 , 4.0 , 5.0 , 11.0 , 8.0 , 9.0 , 5.0 , 7.0 , 7.0 , 2.0 , 8.0 , 7.0 , 4.0 , 5.0 , 10.0 , 12.0 , 8.0 , 10.0 , 6.0 , 10.0 , 6.0 , 10.0 , 12.0 , 11.0 , 14.0 , 12.0 , 23.0 , 22.0 , 29.0 , 54.0 , 108.0 , 298.0 , 1339.0 , 6929.0 , 7292.0 , 6599.0 , 4152.0 , 2058.0 , 1105.0 , 647.0 , 362.0 , 239.0 , 137.0 , 125.0 , 75.0 , 74.0 , 63.0 , 57.0 , 66.0 , 39.0 , 56.0 , 46.0 , 45.0 , 48.0 , 33.0 , 30.0 , 27.0 , 23.0 , 23.0 , 15.0 , 29.0 , 16.0 , 8.0 , 24.0 , 19.0 , 21.0 , 20.0 , 18.0 , 16.0 , 11.0 , 22.0 , 14.0 , 11.0 , 26.0 , 17.0 , 17.0 , 15.0 , 15.0 , 15.0 , 16.0 , 11.0 , 14.0 , 21.0 , 15.0 , 11.0 , 7.0 , 13.0 , 19.0 , 9.0 , 10.0 , 9.0 , 12.0 , 10.0 , 16.0 , 15.0 , 14.0 , 11.0 , 7.0 , 18.0 , 9.0 , 20.0 , 9.0 , 12.0 , 17.0 , 16.0 , 10.0 , 11.0 , 10.0 , 12.0};

  std::vector<float> bins= {-29.8 , -29.4 , -29.0 , -28.6 , -28.2 , -27.8 , -27.4 , -27.0 , -26.6 , -26.2 , -25.8 , -25.4 , -25.0 , -24.6 , -24.2 , -23.8 , -23.4 , -23.0 , -22.6 , -22.2 , -21.8 , -21.4 , -21.0 , -20.6 , -20.2 , -19.8 , -19.4 , -19.0 , -18.6 , -18.2 , -17.8 , -17.4 , -17.0 , -16.6 , -16.2 , -15.8 , -15.4 , -15.0 , -14.6 , -14.2 , -13.8 , -13.4 , -13.0 , -12.6 , -12.2 , -11.8 , -11.4 , -11.0 , -10.6 , -10.2 , -9.8 , -9.4 , -9.0 , -8.6 , -8.2 , -7.8 , -7.4 , -7.0 , -6.6 , -6.2 , -5.8 , -5.4 , -5.0 , -4.6 , -4.2 , -3.8 , -3.4 , -3.0 , -2.6 , -2.2 , -1.8 , -1.4 , -1.0 , -0.6 , -0.2 , 0.2 , 0.6 , 1.0 , 1.4 , 1.8 , 2.2 , 2.6 , 3.0 , 3.4 , 3.8 , 4.2 , 4.6 , 5.0 , 5.4 , 5.8 , 6.2 , 6.6 , 7.0 , 7.4 , 7.8 , 8.2 , 8.6 , 9.0 , 9.4 , 9.8 , 10.2 , 10.6 , 11.0 , 11.4 , 11.8 , 12.2 , 12.6 , 13.0 , 13.4 , 13.8 , 14.2 , 14.6 , 15.0 , 15.4 , 15.8 , 16.2 , 16.6 , 17.0 , 17.4 , 17.8 , 18.2 , 18.6 , 19.0 , 19.4 , 19.8 , 20.2 , 20.6 , 21.0 , 21.4 , 21.8 , 22.2 , 22.6 , 23.0 , 23.4 , 23.8 , 24.2 , 24.6 , 25.0 , 25.4 , 25.8 , 26.2 , 26.6 , 27.0 , 27.4 , 27.8 , 28.2 , 28.6 , 29.0 , 29.4 , 29.8};

  _hist_x = new TH1D("XRes","XRes",100,-20,20);
  _hist_y = new TH1D("YRes","YRes",100,-20,20);
  _hist_z = new TH1D("ZRes","ZRes",60,-20,20);

  for(int i = 0; i < 60; i++){
    _hist_x->SetBinContent(i,x_res[i]); 
    _hist_y->SetBinContent(i,y_res[i]); 
    _hist_z->SetBinContent(i,z_res[i]); 
  }

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
    //std::random_device rd;
    //std::mt19937 e2(rd());
    //std::normal_distribution<double> x_vtx(0,3.75); //2.77); commented out are the pi0 resolutions
    //std::normal_distribution<double> y_vtx(0,3.77); //2.38); comm
    //std::normal_distribution<double> z_vtx(0,3.71); //1.97);
    //double xnew = x_vtx(e2); //generator);
    //double ynew = y_vtx(e2); //generator);
    //double znew = z_vtx(e2); //generator);

    double xnew = _hist_x->GetRandom();
    double ynew = _hist_y->GetRandom();
    double znew = _hist_z->GetRandom();

    std::cout<<"Vertex: "<<xnew<<", "<<ynew<<", "<<znew<<std::endl ;
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
