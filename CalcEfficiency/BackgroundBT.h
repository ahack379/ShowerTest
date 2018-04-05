/**
 * \file BackgroundBT.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class BackgroundBT
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_BACKGROUNDBT_H
#define LARLITE_BACKGROUNDBT_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"


namespace larlite {
  /**
     \class BackgroundBT
     User custom analysis class made by SHELL_USER_NAME
   */
  class BackgroundBT : public ana_base{
  
  public:

    /// Default constructor
    BackgroundBT(){ _name="BackgroundBT"; _fout=0; _tree=0; _shower_tree = 0; _univ=0; _get_genie_info=false; _eventweight_producer = "";}

    /// Default destructor
    virtual ~BackgroundBT(){}

    /** IMPLEMENT in BackgroundBT.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in BackgroundBT.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in BackgroundBT.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void UseMCSample ( bool useit=false ) { _mc_sample = useit; }

    void GetPi0Info  ( bool getit=false ) { _get_pi0_info = getit; }

    void GetSingleShowerInfo  ( bool getit=false ) { _get_single_shower_info = getit; }

    void clear() ;

    // The rest is for uncertainty re-weighting
    // Tell the module to calculate uncertainties 
    void GetUncertaintyInfo  ( bool getit=false ) { _get_genie_info = getit; }

    // Pick the uncertainty to vary (genie or flux)
    void SetEWProducer( std::string producer ) { _eventweight_producer = producer ; }

    // This is for the flux uncertainty portion of the analysis
    void SetNominal ( float nom ) { _N = nom ; }
    void SetNominalXsec ( float nom ) { _N_xsec = nom ; }

  protected:

  std::vector<int> _event_list ;

  int _n_noise;

  int _n_cosmic ;
  int _n_cc1pi0 ; // This is the signal
  int _n_nc1pi0 ; 

  int _n_cc1pi0_outFV ; // This is the signal
  int _n_multpi0 ;
  int _n_cccex ;
  int _n_nccex ;
  int _n_nue ;
  int _n_antimu ;
  int _n_ncother;
  int _n_ccother ;

  int _n_gamma ;
  int _n_kaondecay;
  int _n_other ;

  // One entry per shower
  TTree * _shower_tree ;
  float _shr_startx;
  float _shr_starty;
  float _shr_startz;
  float _shr_startw;
  float _shr_startt;
  float _shr_dirx;
  float _shr_diry;
  float _shr_dirz;
  float _shr_energy;
  float _shr_oangle;
  float _shr_dedx;
  float _shr_vtx_dist;
  float _shr_trk_delta_theta ;
  float _shr_trk_delta_phi ;

  float _shr_purity ;
  float _shr_complete ;
  float _shr_origin ;
  float _shr_type;
  bool _shr_from_pi0 ;

  TTree * _tree ;
  int _event ;
  int _bkgd_id ;
  int _nu_mode ;
  int _nshrs;
  float _vtx_x ;
  float _vtx_y ;
  float _vtx_z ;
  float _mc_vtx_x ;
  float _mc_vtx_y ;
  float _mc_vtx_z ;

  std::map<int,int> _mc_hit_map ;

  float _pi0_mass ;
  float _pi0_mom ;
  float _pi0_oangle ;
  float _pi0_true_oangle;   // is this mccluster from a pi0? yes(1) no(0)
  float _pi0_IP;
  float _pi0_low_shrE;
  float _pi0_high_shrE;
  float _pi0_low_radL;
  float _pi0_high_radL;
  float _pi0_low_IP_w_vtx;
  float _pi0_high_IP_w_vtx;

  float _pi0_low_purity;
  float _pi0_high_purity;
  float _pi0_low_complete;
  float _pi0_high_complete;
  float _pi0_low_cw_purity;
  float _pi0_high_cw_purity;
  float _pi0_low_cw_complete;
  float _pi0_high_cw_complete;
  float _pi0_low_true_gammaE;
  float _pi0_high_true_gammaE;
  float _pi0_low_true_detProf_gammaE;
  float _pi0_high_true_detProf_gammaE;
  float _pi0_low_reco_gammaE;
  float _pi0_high_reco_gammaE;
  // Also add some other truth info from the mccluster
  float _pi0_low_origin ; // is the corresponding mccluster due to nu(1), cosmic(2), noise(3)
  float _pi0_low_type ;   // is this mccluster due to track(0) or shower(1)
  bool _pi0_low_from_pi0;   // is this mccluster from a pi0? yes(1) no(0)
  float _pi0_high_origin ; // is this mccluster due to nu(1), cosmic(2), noise(3)
  float _pi0_high_type ;   // is this mccluster due to track(0) or shower(1)
  bool _pi0_high_from_pi0;   // is this mccluster from a pi0? yes(1) no(0)
  
  float _pi0_low_dist_to_nearest_trk ;
  float _pi0_high_dist_to_nearest_trk ;

  float _pi0_low_st_x ;
  float _pi0_low_st_y ;
  float _pi0_low_st_z ;
  float _pi0_high_st_x ;
  float _pi0_high_st_y ;
  float _pi0_high_st_z ;
  float _pi0_low_true_st_x ;
  float _pi0_low_true_st_y ;
  float _pi0_low_true_st_z ;
  float _pi0_high_true_st_x ;
  float _pi0_high_true_st_y ;
  float _pi0_high_true_st_z ;


  //////////////////////////////////////
  float _gamma_E;
  float _gamma_RL ;
  float _gamma_IP_w_vtx;
  float _gamma_purity; 
  float _gamma_complete; 
  float _gamma_cw_purity; 
  float _gamma_cw_complete; 
  float _gamma_trueE;
  float _gamma_trueE_detProf;
  // Also add some other truth info from the mccluster
  float _gamma_origin ; // is this mccluster due to nu(1), cosmic(2), noise(-1)
  float _gamma_type ;   // is this mccluster due to track(0) or shower(1)
  bool _gamma_from_pi0;   // is this mccluster from a pi0? yes(1) no(0)


  float _mu_angle ;
  float _mu_len ;
  float _mu_startx ;
  float _mu_starty ;
  float _mu_startz ;
  float _mu_endx ;
  float _mu_endy ;
  float _mu_endz ;
  float _mu_phi ;
  float _mu_mom;
  float _mult;
  float _mu_purity ;
  float _mu_complete ;
  float _mu_cw_purity ;
  float _mu_cw_complete ;
  // Also add some other truth info from the mccluster
  float _mu_origin ; // is this mccluster due to nu(1), cosmic(2), noise(3)
  float _mu_type ;   // is this mccluster due to track(0) or shower(1)

  int _n_shr_pi0 ;
  int _n_shr_nushr ;
  int _n_shr_nutrk ;
  int _n_shr_cosmic ;
  int _n_shr_noise ;

  // Flashes
  float _flash_time ;
  int _flash_pe ;
  float _flash_y_center ;
  float _flash_z_center ;
  float _flash_y_width;
  float _flash_z_width ;

  ::geoalgo::GeoAlgo _geoAlgo ;

  bool _mc_sample ;
  bool _get_pi0_info ;
  bool _get_single_shower_info ;
  bool _get_genie_info ;

  std::vector<std::string> _bkgd_v ;

  larutil::SpaceChargeMicroBooNE *_SCE;
  double _time2cm;

  // Post technote version v0.9 additions
  // "tree" additions
  int _n_track_hits_0 ;
  int _n_track_hits_1 ;
  int _n_track_hits_2 ;
  int _n_shower_hits_0 ;
  int _n_shower_hits_1 ;
  int _n_shower_hits_2 ;

  // "showertree" additions
  float _shr_ip ;
  float _shr_rl;

  // Check osc group question
  int _1gamma ;

  // genie addition
  std::vector<float> _sel_evts_m1;
  std::vector<float> _sel_evts_p1;

  std::vector<std::string> _genie_label_v ;
  std::string _eventweight_producer ;

  // FLUX variables
  std::vector<std::string> _unisim_label_v ;
  bool _signal ;

  std::vector<std::vector<float>>  _s_weights_by_universe ;
  std::vector<std::vector<float>>  _b_weights_by_universe ;
  std::vector<std::vector<float>>  _e_weights_by_universe ;
  std::vector<std::vector<float>>  _t_weights_by_universe ;
  std::vector<std::vector<float>>  _flux_by_universe ; 

  TTree * _univ;
  std::vector<std::vector<float>> _xsec_v ;
  std::vector<std::vector<float>> _perc_v;

  std::map<std::string,int> _label_map ;

  float _N ;
  float _N_xsec ;

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
