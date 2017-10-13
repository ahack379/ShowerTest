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

namespace larlite {
  /**
     \class BackgroundBT
     User custom analysis class made by SHELL_USER_NAME
   */
  class BackgroundBT : public ana_base{
  
  public:

    /// Default constructor
    BackgroundBT(){ _name="BackgroundBT"; _fout=0; _tree=0; }

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

  //////////////////////////////////////
  float _gamma_E;
  float _gamma_RL ;
  float _gamma_purity; 
  float _gamma_complete; 
  float _gamma_cw_purity; 
  float _gamma_cw_complete; 

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

  ::geoalgo::GeoAlgo _geoAlgo ;

  bool _mc_sample ;
  bool _get_pi0_info ;
  bool _get_single_shower_info ;

  std::vector<std::string> _bkgd_v ;
    
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
