/**
 * \file LowLevel.h
 *
 * \ingroup PlotDump
 * 
 * \brief Class def header for a class LowLevel
 *
 * @author ah673
 */

/** \addtogroup PlotDump

    @{*/

#ifndef LARLITE_LOWLEVEL_H
#define LARLITE_LOWLEVEL_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class LowLevel
     User custom analysis class made by SHELL_USER_NAME
   */
  class LowLevel : public ana_base{
  
  public:

    /// Default constructor
    LowLevel(){ _name="LowLevel"; _fout=0; _low_level_tree = 0; _hit_tree = 0; _shower_tree = 0;}

    /// Default destructor
    virtual ~LowLevel(){}

    /** IMPLEMENT in LowLevel.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in LowLevel.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in LowLevel.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void clear();

    void UseMCSample( bool usemc=true ){ _mc_sample = usemc ; }

  protected:
    
    // One entry per event
    TTree * _low_level_tree ;
    int _nhits0;
    int _nhits1;
    int _nhits2;

    float _startx;
    float _starty;
    float _startz;
    float _endx;
    float _endy;
    float _endz;
    float _length;
    float _theta;
    float _phi;
    float _mult;
    
    float _vtxx;
    float _vtxy;
    float _vtxz;
    float _vtxw;
    float _vtxt;
    float _vtx_trk_dist ;
    float _vtx_mc_reco_dist;
    float _mc_vtxx;
    float _mc_vtxy;
    float _mc_vtxz;
    float _mc_vtxw;
    float _mc_vtxt;

    int _nshrs; 
    
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

    // One entry per hit 
    TTree * _hit_tree ;
    int _plane ;
    float _charge;
    float _wire ;
    float _time_peak;
    float _time_width;
    float _gof ;
    
    // One entry per shower
    TTree * _shower_tree ;

    int _entry ;
    bool _mc_sample;
    
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
