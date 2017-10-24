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
    BackgroundBT(){ _name="BackgroundBT"; _fout=0; _tree=0; _shower_tree = 0; }

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
  float _pi0_IP;
  float _pi0_low_shrE;
  float _pi0_high_shrE;
  float _pi0_low_radL;
  float _pi0_high_radL;

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

  //////////////////////////////////////
  float _gamma_E;
  float _gamma_RL ;
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

  ::geoalgo::GeoAlgo _geoAlgo ;

  bool _mc_sample ;
  bool _get_pi0_info ;
  bool _get_single_shower_info ;

  std::vector<std::string> _bkgd_v ;

  larutil::SpaceChargeMicroBooNE *_SCE;
  double _time2cm;

    
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
