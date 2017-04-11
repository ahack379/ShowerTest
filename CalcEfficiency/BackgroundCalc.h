/**
 * \file BackgroundCalc.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class BackgroundCalc
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_BACKGROUNDCALC_H
#define LARLITE_BACKGROUNDCALC_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class BackgroundCalc
     User custom analysis class made by SHELL_USER_NAME
   */
  class BackgroundCalc : public ana_base{
  
  public:

    /// Default constructor
    BackgroundCalc(){ _name="BackgroundCalc"; _fout=0; _tree=0; }

    /// Default destructor
    virtual ~BackgroundCalc(){}

    /** IMPLEMENT in BackgroundCalc.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in BackgroundCalc.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in BackgroundCalc.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  std::vector<int> _event_list ;
  std::vector<int> _ccother_list ;
  std::vector<int> _pi0_list ;

  int _n_noise ;
  int _n_cosmic ;
  int _n_nue ;
  int _n_antinumu ;
  int _n_nc;
  int _multpi0 ;
  int _ccpi0_outfv;
  int _tot_ccpi0; // This is the signal
  int _n_gammas;
  int _n_ccother ;

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

  ::geoalgo::GeoAlgo _geoAlgo ;
    
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
