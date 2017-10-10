/**
 * \file BackgroundObnoxious.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class BackgroundObnoxious
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_BACKGROUNDOBNOXIOUS_H
#define LARLITE_BACKGROUNDOBNOXIOUS_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class BackgroundObnoxious
     User custom analysis class made by SHELL_USER_NAME
   */
  class BackgroundObnoxious : public ana_base{
  
  public:

    /// Default constructor
    BackgroundObnoxious(){ _name="BackgroundObnoxious"; _fout=0; _tree=0; }

    /// Default destructor
    virtual ~BackgroundObnoxious(){}

    /** IMPLEMENT in BackgroundObnoxious.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in BackgroundObnoxious.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in BackgroundObnoxious.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void UseMCSample ( bool useit=false ) { _mc_sample = useit; }

    void GetPi0Info  ( bool getit=false ) { _get_pi0_info = getit; }

    void GetSingleShowerInfo  ( bool getit=false ) { _get_single_shower_info = getit; }

    void clear() ;

  protected:

  std::vector<int> _event_list ;
  std::vector<int> _vtx_list ;

  int _n_other ;
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


  TTree * _tree ;
  int _event ;
  int _bkgd_id ;
  int _nu_mode ;
  int _nshrs;
  float _vtx_x ;
  float _vtx_y ;
  float _vtx_z ;

  float _pi0_mass ;
  float _pi0_mom ;
  float _pi0_oangle ;
  float _pi0_IP;
  float _pi0_low_shrE;
  float _pi0_high_shrE;
  float _pi0_low_radL;
  float _pi0_high_radL;

  float _gamma_E;
  float _gamma_RL ;

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
