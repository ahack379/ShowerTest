/**
 * \file ShowerQuality.h
 *
 * \ingroup Signal_Quality
 * 
 * \brief Class def header for a class ShowerQuality
 *
 * @author ah673
 */

/** \addtogroup Signal_Quality

    @{*/

#ifndef LARLITE_SHOWERQUALITY_H
#define LARLITE_SHOWERQUALITY_H

#include "Analysis/ana_base.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class ShowerQuality
     User custom analysis class made by SHELL_USER_NAME
   */
  class ShowerQuality : public ana_base{
  
  public:

    /// Default constructor
    ShowerQuality(){ _name="ShowerQuality"; _fout=0; _vtx_tree=0; _offset = 0.7; }

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

    void SetXOffset(double d) { _offset = d; }

  protected:

  std::multimap<float,float> _map;
  std::vector<int> _event_list ;
  int _event ; 

  TTree * _vtx_tree;

  float _vtx_diff ;
  float _mc_vtx_x ;
  float _mc_vtx_y ;
  float _mc_vtx_z ;
  float _reco_vtx_x ;
  float _reco_vtx_y ;
  float _reco_vtx_z ;

  // Track variables
  float _trk_st_diff ;
  float _reco_trk_st_x ;
  float _reco_trk_st_y ;
  float _reco_trk_st_z ;
  float _mc_trk_dir_x ;
  float _mc_trk_dir_y ;
  float _mc_trk_dir_z ;
  float _reco_trk_dir_x ;
  float _reco_trk_dir_y ;
  float _reco_trk_dir_z ;
  float _trk_dot ;

  float _reco_shr1_e;
  float _mc_shr1_e;
  float _shr1_dot; 
  float _reco_shr1_st_x;
  float _reco_shr1_st_y;
  float _reco_shr1_st_z;
  double _mc_shr1_st_x;
  double _mc_shr1_st_y;
  double _mc_shr1_st_z;
  float _shr1_st_diff;

  float _mc_shr2_e;
  float _reco_shr2_e;
  float _shr2_dot; 
  float _reco_shr2_st_x;
  float _reco_shr2_st_y;
  float _reco_shr2_st_z;
  double _mc_shr2_st_x;
  double _mc_shr2_st_y;
  double _mc_shr2_st_z;
  float _shr2_st_diff;

  double _offset;
  double _time2cm;
  larutil::SpaceChargeMicroBooNE *_SCE;

  geoalgo::GeoAlgo _geoAlgo ;

  int _bad_events ;
  int _pi0s; 
  int _fv ;

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
