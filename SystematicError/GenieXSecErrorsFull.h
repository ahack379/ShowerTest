/**
 * \file GenieXSecErrorsFull.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class GenieXSecErrorsFull
 *
 * @author ariana hackenburg 
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_GENIEXSECERRORSFULL_H
#define LARLITE_GENIEXSECERRORSFULL_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"
#include <fstream>

namespace larlite {
  /**
     \class GenieXSecErrorsFull
     User custom analysis class made by SHELL_USER_NAME
   */
  class GenieXSecErrorsFull : public ana_base{
  
  public:

    /// Default constructor
    GenieXSecErrorsFull(); //{ _name="GenieXSecErrorsFull"; _fout=0;}

    /// Default destructor
    virtual ~GenieXSecErrorsFull(){}

    /** IMPLEMENT in GenieXSecErrorsFull.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GenieXSecErrorsFull.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GenieXSecErrorsFull.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    std::ifstream _file; //("list.txt",std::ifstream::in);

    float _tot_pot ;
    const larutil::Geometry * fGeometry ;

    std::multimap<float,float> _map ;
    
    /////////// My extra variables
    TTree * _tree ;
    float _xsec_mom_truth ;
    float _xsec_theta_truth ;
    std::vector<float> _weight_v ;

    TTree * _final_tree ;
    float _all_evts_nominal;
    std::vector<float> _all_evts_m1 ;
    std::vector<float> _all_evts_p1 ;
    std::vector<std::string> _genie_label_v ;


    //std::vector<int> _tagged_v ;
    //std::vector<int> _bkgd_v ;
    //std::vector<int> _eff_v ;
    //std::vector<int> _xsec_v ;

    //std::vector<float> _genie_evwgt_v;

    //std::vector<float> _xsec_mutheta_v ;
    //std::vector<float> _xsec_pi0theta_v ;
    //std::vector<float> _xsec_pi0mom_v ;

 
    //double all_evts_nominal = 0;  // all numu cc events (w/o selection) -> nominal values of the parameters
    //vector<double> all_evts_p1;   // all numu cc events (w/o selection) -> + 1 sigma parameters
    //vector<double> all_evts_m1;   // all numu cc events (w/o selection) -> - 1 sigma parameters
    //
    //double sel_evts_nominal = 0;  // selected numu cc events (w/o selection) -> nominal values of the parameters
    //vector<double> sel_evts_p1;   // selected numu cc events (w/o selection) -> + 1 sigma parameters
    //vector<double> sel_evts_m1;   // selected numu cc events (w/o selection) -> - 1 sigma parameters
    
    //TH1F *xsec_mom_truth = new TH1F("xsec_mom_truth", "", 10, 0, 2);
    //TH1F *xsec_mom_data = new TH1F("xsec_mom_data", "", 10, 0, 2);
    //TH1F *xsec_mom_bg = new TH1F("xsec_mom_bg", "", 10, 0, 2);
    //TH1F *xsec_mom_eff = new TH1F("xsec_mom_eff", "", 10, 0, 2);
    //
    //TH1F *xsec_mom_reco_truth = new TH1F("xsec_mom_reco_truth", "", 10, 0, 2);
    //TH1F *xsec_mom_reco_data = new TH1F("xsec_mom_reco_data", "", 10, 0, 2);
    //TH1F *xsec_mom_reco_bg = new TH1F("xsec_mom_reco_bg", "", 10, 0, 2);
    //TH1F *xsec_mom_reco_eff = new TH1F("xsec_mom_reco_eff", "", 10, 0, 2);
    //
    //TH1F *xsec_theta_truth = new TH1F("xsec_theta_truth", "", 10, -1, 1);
    //TH1F *xsec_theta_data = new TH1F("xsec_theta_data", "", 10, -1, 1);
    //TH1F *xsec_theta_bg = new TH1F("xsec_theta_bg", "", 10, -1, 1);
    //TH1F *xsec_theta_eff = new TH1F("xsec_theta_eff", "", 10, -1, 1);
    //
    //TH1F *xsec_theta_reco_truth = new TH1F("xsec_theta_reco_truth", "", 10, -1, 1);
    //TH1F *xsec_theta_reco_data = new TH1F("xsec_theta_reco_data", "", 10, -1, 1);
    //TH1F *xsec_theta_reco_bg = new TH1F("xsec_theta_reco_bg", "", 10, -1, 1);
    //TH1F *xsec_theta_reco_eff = new TH1F("xsec_theta_reco_eff", "", 10, -1, 1);
    //
    //vector<TH1D*> xsec_mom_truth_p1;
    //vector<TH1D*> xsec_mom_truth_m1;
    //vector<TH1D*> xsec_theta_truth_p1;
    //vector<TH1D*> xsec_theta_truth_m1;
    //vector<TH1D*> xsec_mom_data_p1;
    //vector<TH1D*> xsec_mom_data_m1;
    //vector<TH1D*> xsec_mom_reco_data_p1;
    //vector<TH1D*> xsec_mom_reco_data_m1;
    //vector<TH1D*> xsec_theta_data_p1;
    //vector<TH1D*> xsec_theta_data_m1;
    //vector<TH1D*> xsec_theta_reco_data_p1;
    //vector<TH1D*> xsec_theta_reco_data_m1;
    //
    //vector<TH1D*> xsec_mom_eff_p1;
    //vector<TH1D*> xsec_mom_eff_m1;
    //vector<TH1D*> xsec_mom_reco_eff_p1;
    //vector<TH1D*> xsec_mom_reco_eff_m1;
    //vector<TH1D*> xsec_theta_eff_p1;
    //vector<TH1D*> xsec_theta_eff_m1;
    //vector<TH1D*> xsec_theta_reco_eff_p1;
    //vector<TH1D*> xsec_theta_reco_eff_m1;
    //
    //
    //TH1F *pmu_numu_cc_reco_histo = new TH1F("pmu_numu_cc_reco_histo", "", 10, 0, 2);
    //TH1F *costhetamu_numu_cc_reco_histo = new TH1F("costhetamu_numu_cc_reco_histo", "", 10, -1, 1);
    //TH1F *pmu_nc_reco_histo = new TH1F("pmu_nc_reco_histo", "", 10, 0, 2);
    //TH1F *costhetamu_nc_reco_histo = new TH1F("costhetamu_nc_reco_histo", "", 10, -1, 1);
    //
    //
    //vector<TH1D*> pmu_numu_cc_reco_histo_p1;
    //vector<TH1D*> pmu_numu_cc_reco_histo_m1;
    //vector<TH1D*> costhetamu_numu_cc_reco_histo_p1;
    //vector<TH1D*> costhetamu_numu_cc_reco_histo_m1;
    //vector<TH1D*> costhetamu_nue_cc_reco_histo_p1;
    //vector<TH1D*> costhetamu_nue_cc_reco_histo_m1;
    //vector<TH1D*> pmu_nc_reco_histo_p1;
    //vector<TH1D*> pmu_nc_reco_histo_m1;
    //vector<TH1D*> costhetamu_nc_reco_histo_p1;
    //vector<TH1D*> costhetamu_nc_reco_histo_m1;
    //
    //
    //// Cross Section
    //TH1D *        XSec_pmu_nominal;
    //vector<TH1D*> XSec_pmu_p1;
    //vector<TH1D*> XSec_pmu_m1;
    //
    //TH1D * background_pmu_nominal;
    //vector<TH1D*> background_pmu_p1;
    //vector<TH1D*> background_pmu_m1;
    //
    //vector<TH1D*> XSec_pmu_percDiff_p1;
    //vector<TH1D*> XSec_pmu_percDiff_m1;
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
