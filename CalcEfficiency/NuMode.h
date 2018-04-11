/**
 * \file NuMode.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class NuMode
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_NUMODE_H
#define LARLITE_NUMODE_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class NuMode
     User custom analysis class made by SHELL_USER_NAME
   */
  class NuMode : public ana_base{
  
  public:

    /// Default constructor
    NuMode(){ _name="NuMode"; _fout=0; _tree=0; _shower_tree = 0; _get_gt0_shower = false; }

    /// Default destructor
    virtual ~NuMode(){}

    /** IMPLEMENT in NuMode.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in NuMode.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in NuMode.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void UseMCSample ( bool useit=false ) { _mc_sample = useit; }

    void GetPi0Info  ( bool getit=false ) { _get_pi0_info = getit; }

    void GetSingleShowerInfo  ( bool getit=false ) { _get_single_shower_info = getit; }

    void GetGT0Shower( bool getit=false ) { _get_gt0_shower = getit; }

    void clear() ;

  protected:

  std::vector<int> _event_list ;

  int _n_other ;
  int _n_cosmic ;
  int _n_cc1pi0 ; // This is the signal
  int _n_cc0pi0 ;
  int _n_nc1pi0 ; 
  int _n_nc0pi0 ;

  TTree * _tree ;
  int _event ;
  int _bkgd_id ;
  int _nu_mode;
  int _nshrs ;
  float _vtx_x ;
  float _vtx_y ;
  float _vtx_z ;

  float _pi0_mass ;
  float _pi0_mom ;
  float _pi0_oangle ;
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


  ::geoalgo::GeoAlgo _geoAlgo ;

  bool _mc_sample ;
  bool _get_pi0_info ;
  bool _get_single_shower_info ;

  bool _get_gt0_shower ;

  std::vector<std::multimap<float,float>> _map_v;
  
    
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
