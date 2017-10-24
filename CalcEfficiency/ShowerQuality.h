/**
 * \file ShowerQuality.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class ShowerQuality
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_SHOWERQUALITY_H
#define LARLITE_SHOWERQUALITY_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

#include "LArUtil/SpaceChargeMicroBooNE.h"

namespace larlite {
  /**
     \class ShowerQuality
     User custom analysis class made by SHELL_USER_NAME
   */
  class ShowerQuality : public ana_base{
  
  public:

    /// Default constructor
    ShowerQuality(){ _name="ShowerQuality"; _fout=0; _tree=0; }

    /// Default destructor
    virtual ~ShowerQuality(){}

    /** IMPLEMENT in ShowerQuality.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ShowerQuality.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ShowerQuality.cc! 
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

  float _purity;
  float _complete;
  float _cw_purity;
  float _cw_complete;
  float _mc_e;
  float _mc_detProf_e;
  float _mc_st_x ;
  float _mc_st_y ;
  float _mc_st_z ;

  float _origin ; // is the corresponding mccluster due to nu(1), cosmic(2), noise(3)
  float _type ;   // is this mccluster due to track(0) or shower(1)
  bool _from_pi0;   // is this mccluster from a pi0? yes(1) no(0)

  float _reco_e;
  float _st_x ;
  float _st_y ;
  float _st_z ;

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
