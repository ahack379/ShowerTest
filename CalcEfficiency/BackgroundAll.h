/**
 * \file BackgroundAll.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class BackgroundAll
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_BACKGROUNDALL_H
#define LARLITE_BACKGROUNDALL_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class BackgroundAll
     User custom analysis class made by SHELL_USER_NAME
   */
  class BackgroundAll : public ana_base{
  
  public:

    /// Default constructor
    BackgroundAll(){ _name="BackgroundAll"; _fout=0; _tree=0; }

    /// Default destructor
    virtual ~BackgroundAll(){}

    /** IMPLEMENT in BackgroundAll.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in BackgroundAll.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in BackgroundAll.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void UseMCSample ( bool useit=false ) { _mc_sample = useit; }

    void GetPi0Info  ( bool getit=false ) { _get_pi0_info = getit; }

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