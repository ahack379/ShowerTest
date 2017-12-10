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

  std::vector<float> x_res={2.0 , 14.0 , 4.0 , 5.0 , 2.0 , 5.0 , 6.0 , 10.0 , 4.0 , 5.0 , 8.0 , 5.0 , 11.0 , 6.0 , 4.0 , 5.0 , 6.0 , 11.0 , 5.0 , 11.0 , 8.0 , 9.0 , 4.0 , 9.0 , 13.0 , 9.0 , 11.0 , 4.0 , 7.0 , 7.0 , 8.0 , 7.0 , 6.0 , 3.0 , 8.0 , 11.0 , 6.0 , 8.0 , 6.0 , 10.0 , 9.0 , 8.0 , 6.0 , 9.0 , 10.0 , 8.0 , 7.0 , 10.0 , 8.0 , 11.0 , 11.0 , 10.0 , 5.0 , 9.0 , 6.0 , 9.0 , 11.0 , 6.0 , 11.0 , 15.0 , 13.0 , 15.0 , 10.0 , 9.0 , 16.0 , 13.0 , 8.0 , 12.0 , 16.0 , 14.0 , 8.0 , 13.0 , 12.0 , 17.0 , 15.0 , 6.0 , 14.0 , 18.0 , 11.0 , 10.0 , 13.0 , 19.0 , 22.0 , 23.0 , 24.0 , 21.0 , 25.0 , 25.0 , 28.0 , 32.0 , 42.0 , 43.0 , 64.0 , 76.0 , 123.0 , 239.0 , 619.0 , 1540.0 , 3179.0 , 7331.0 , 10827.0 , 4372.0 , 2073.0 , 860.0 , 348.0 , 150.0 , 79.0 , 63.0 , 49.0 , 35.0 , 22.0 , 26.0 , 25.0 , 24.0 , 26.0 , 25.0 , 16.0 , 22.0 , 16.0 , 9.0 , 15.0 , 15.0 , 14.0 , 16.0 , 16.0 , 12.0 , 10.0 , 12.0 , 11.0 , 12.0 , 16.0 , 12.0 , 13.0 , 10.0 , 8.0 , 9.0 , 8.0 , 9.0 , 14.0 , 5.0 , 7.0 , 6.0 , 8.0 , 7.0 , 14.0 , 13.0 , 7.0 , 7.0 , 8.0 , 8.0 , 8.0 , 10.0 , 5.0 , 7.0 , 10.0 , 11.0 , 12.0 , 12.0 , 6.0 , 14.0 , 7.0 , 2.0 , 10.0 , 8.0 , 10.0 , 6.0 , 14.0 , 8.0 , 11.0 , 4.0 , 12.0 , 3.0 , 5.0 , 9.0 , 11.0 , 7.0 , 5.0 , 2.0 , 6.0 , 4.0 , 8.0 , 2.0 , 13.0 , 10.0 , 5.0 , 3.0 , 11.0 , 7.0 , 6.0 , 10.0 , 14.0 , 6.0 , 7.0 , 5.0 , 4.0 , 5.0 , 7.0 , 8.0 , 11.0 , 8.0};

  std::vector<float> y_res = {7.0 , 2.0 , 3.0 , 7.0 , 5.0 , 6.0 , 7.0 , 8.0 , 6.0 , 3.0 , 7.0 , 9.0 , 3.0 , 7.0 , 5.0 , 5.0 , 3.0 , 7.0 , 6.0 , 1.0 , 10.0 , 7.0 , 5.0 , 11.0 , 9.0 , 11.0 , 10.0 , 11.0 , 8.0 , 7.0 , 10.0 , 12.0 , 6.0 , 13.0 , 4.0 , 10.0 , 14.0 , 10.0 , 6.0 , 13.0 , 9.0 , 11.0 , 18.0 , 10.0 , 11.0 , 14.0 , 9.0 , 11.0 , 8.0 , 9.0 , 11.0 , 17.0 , 12.0 , 9.0 , 15.0 , 9.0 , 11.0 , 18.0 , 9.0 , 17.0 , 8.0 , 14.0 , 21.0 , 10.0 , 19.0 , 18.0 , 10.0 , 19.0 , 14.0 , 16.0 , 18.0 , 16.0 , 19.0 , 18.0 , 16.0 , 17.0 , 34.0 , 23.0 , 19.0 , 23.0 , 24.0 , 21.0 , 19.0 , 25.0 , 30.0 , 26.0 , 26.0 , 46.0 , 60.0 , 33.0 , 61.0 , 70.0 , 104.0 , 166.0 , 273.0 , 444.0 , 805.0 , 1624.0 , 3678.0 , 9382.0 , 9477.0 , 3293.0 , 1301.0 , 669.0 , 342.0 , 210.0 , 131.0 , 101.0 , 68.0 , 56.0 , 45.0 , 33.0 , 29.0 , 24.0 , 35.0 , 41.0 , 28.0 , 27.0 , 12.0 , 22.0 , 30.0 , 17.0 , 17.0 , 25.0 , 18.0 , 15.0 , 16.0 , 16.0 , 12.0 , 20.0 , 12.0 , 13.0 , 18.0 , 9.0 , 13.0 , 13.0 , 13.0 , 8.0 , 9.0 , 13.0 , 11.0 , 16.0 , 15.0 , 8.0 , 9.0 , 10.0 , 9.0 , 7.0 , 10.0 , 14.0 , 12.0 , 17.0 , 10.0 , 12.0 , 4.0 , 8.0 , 11.0 , 5.0 , 15.0 , 12.0 , 8.0 , 5.0 , 10.0 , 6.0 , 7.0 , 13.0 , 8.0 , 11.0 , 7.0 , 9.0 , 4.0 , 6.0 , 7.0 , 6.0 , 4.0 , 7.0 , 6.0 , 8.0 , 8.0 , 7.0 , 12.0 , 9.0 , 4.0 , 4.0 , 6.0 , 11.0 , 4.0 , 10.0 , 6.0 , 4.0 , 4.0 , 5.0 , 10.0 , 2.0 , 3.0 , 7.0 , 5.0 , 5.0 , 7.0 , 7.0};

  std::vector<float> z_res = {2.0 , 10.0 , 3.0 , 2.0 , 3.0 , 6.0 , 3.0 , 9.0 , 6.0 , 5.0 , 0.0 , 4.0 , 1.0 , 4.0 , 2.0 , 4.0 , 1.0 , 4.0 , 8.0 , 6.0 , 4.0 , 5.0 , 4.0 , 6.0 , 5.0 , 3.0 , 5.0 , 6.0 , 4.0 , 6.0 , 1.0 , 2.0 , 5.0 , 5.0 , 2.0 , 3.0 , 4.0 , 8.0 , 3.0 , 10.0 , 6.0 , 4.0 , 2.0 , 5.0 , 2.0 , 6.0 , 2.0 , 4.0 , 2.0 , 5.0 , 4.0 , 6.0 , 8.0 , 4.0 , 6.0 , 1.0 , 4.0 , 2.0 , 4.0 , 10.0 , 6.0 , 6.0 , 7.0 , 3.0 , 7.0 , 4.0 , 3.0 , 2.0 , 6.0 , 5.0 , 6.0 , 2.0 , 4.0 , 6.0 , 9.0 , 8.0 , 7.0 , 5.0 , 8.0 , 4.0 , 10.0 , 1.0 , 6.0 , 9.0 , 9.0 , 9.0 , 8.0 , 11.0 , 10.0 , 14.0 , 16.0 , 17.0 , 23.0 , 25.0 , 48.0 , 95.0 , 179.0 , 564.0 , 1751.0 , 6072.0 , 5724.0 , 4928.0 , 4543.0 , 2848.0 , 1647.0 , 1059.0 , 649.0 , 455.0 , 297.0 , 197.0 , 146.0 , 98.0 , 96.0 , 67.0 , 54.0 , 57.0 , 52.0 , 43.0 , 40.0 , 51.0 , 30.0 , 35.0 , 41.0 , 35.0 , 32.0 , 40.0 , 32.0 , 22.0 , 22.0 , 25.0 , 16.0 , 17.0 , 17.0 , 13.0 , 14.0 , 23.0 , 12.0 , 9.0 , 9.0 , 18.0 , 11.0 , 22.0 , 12.0 , 15.0 , 12.0 , 15.0 , 10.0 , 8.0 , 15.0 , 15.0 , 9.0 , 8.0 , 16.0 , 22.0 , 9.0 , 13.0 , 10.0 , 12.0 , 10.0 , 13.0 , 12.0 , 9.0 , 8.0 , 12.0 , 20.0 , 11.0 , 6.0 , 10.0 , 4.0 , 8.0 , 12.0 , 15.0 , 5.0 , 8.0 , 14.0 , 1.0 , 10.0 , 7.0 , 8.0 , 13.0 , 11.0 , 8.0 , 14.0 , 7.0 , 5.0 , 12.0 , 9.0 , 8.0 , 14.0 , 7.0 , 11.0 , 9.0 , 14.0 , 14.0 , 6.0 , 9.0 , 8.0 , 7.0 , 6.0 , 12.0};

std::vector<float>bins={-29.85 , -29.55 , -29.25 , -28.95 , -28.65 , -28.35 , -28.05 , -27.75 , -27.45 , -27.15 , -26.85 , -26.55 , -26.25 , -25.95 , -25.65 , -25.35 , -25.05 , -24.75 , -24.45 , -24.15 , -23.85 , -23.55 , -23.25 , -22.95 , -22.65 , -22.35 , -22.05 , -21.75 , -21.45 , -21.15 , -20.85 , -20.55 , -20.25 , -19.95 , -19.65 , -19.35 , -19.05 , -18.75 , -18.45 , -18.15 , -17.85 , -17.55 , -17.25 , -16.95 , -16.65 , -16.35 , -16.05 , -15.75 , -15.45 , -15.15 , -14.85 , -14.55 , -14.25 , -13.95 , -13.65 , -13.35 , -13.05 , -12.75 , -12.45 , -12.15 , -11.85 , -11.55 , -11.25 , -10.95 , -10.65 , -10.35 , -10.05 , -9.75 , -9.45 , -9.15 , -8.85 , -8.55 , -8.25 , -7.95 , -7.65 , -7.35 , -7.05 , -6.75 , -6.45 , -6.15 , -5.85 , -5.55 , -5.25 , -4.95 , -4.65 , -4.35 , -4.05 , -3.75 , -3.45 , -3.15 , -2.85 , -2.55 , -2.25 , -1.95 , -1.65 , -1.35 , -1.05 , -0.75 , -0.45 , -0.15 , 0.15 , 0.45 , 0.75 , 1.05 , 1.35 , 1.65 , 1.95 , 2.25 , 2.55 , 2.85 , 3.15 , 3.45 , 3.75 , 4.05 , 4.35 , 4.65 , 4.95 , 5.25 , 5.55 , 5.85 , 6.15 , 6.45 , 6.75 , 7.05 , 7.35 , 7.65 , 7.95 , 8.25 , 8.55 , 8.85 , 9.15 , 9.45 , 9.75 , 10.05 , 10.35 , 10.65 , 10.95 , 11.25 , 11.55 , 11.85 , 12.15 , 12.45 , 12.75 , 13.05 , 13.35 , 13.65 , 13.95 , 14.25 , 14.55 , 14.85 , 15.15 , 15.45 , 15.75 , 16.05 , 16.35 , 16.65 , 16.95 , 17.25 , 17.55 , 17.85 , 18.15 , 18.45 , 18.75 , 19.05 , 19.35 , 19.65 , 19.95 , 20.25 , 20.55 , 20.85 , 21.15 , 21.45 , 21.75 , 22.05 , 22.35 , 22.65 , 22.95 , 23.25 , 23.55 , 23.85 , 24.15 , 24.45 , 24.75 , 25.05 , 25.35 , 25.65 , 25.95 , 26.25 , 26.55 , 26.85 , 27.15 , 27.45 , 27.75 , 28.05 , 28.35 , 28.65 , 28.95 , 29.25 , 29.55 , 29.85};

std::cout<<"bins size: "<<bins.size()<<", "<< z_res.size()<<", "<<y_res.size()<<std::endl ;

  _hist_x = new TH1D("XRes","XRes",200,-30,30);
  _hist_y = new TH1D("YRes","YRes",200,-30,30);
  _hist_z = new TH1D("ZRes","ZRes",200,-30,30);

  for(int i = 0; i < 150; i++){
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
