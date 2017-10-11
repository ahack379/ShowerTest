/**
 * \file DataJustifyPi0Cuts.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class DataJustifyPi0Cuts
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_DATAJUSTIFYPI0CUTS_H
#define LARLITE_DATAJUSTIFYPI0CUTS_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/GeometryHelper.h"

namespace larlite {
  /**
     \class DataJustifyPi0Cuts
     User custom analysis class made by SHELL_USER_NAME
   */
  class DataJustifyPi0Cuts : public ana_base{
  
  public:

    /// Default constructor
    DataJustifyPi0Cuts(){ _name="DataJustifyPi0Cuts"; _fout=0; _gamma_tree=0; _tree=0; _one_gamma_tree = 0; }

    /// Default destructor
    virtual ~DataJustifyPi0Cuts(){}

    /** IMPLEMENT in DataJustifyPi0Cuts.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in DataJustifyPi0Cuts.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in DataJustifyPi0Cuts.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void clear();

  protected:

    TTree * _gamma_tree ;
    int _event;
    int _event_type; // signal 0, bkgd with pi0 1, bkgd with no pi0 2
    float _pi0_mass;
    float _pi0_mom;
    float _gamma_oangle;
    float _gamma_low_E;
    float _gamma_high_E;
    float _gamma_low_RL;
    float _gamma_high_RL;
    float _gamma_IP;
    float _gamma1_vtx_IP ;
    float _gamma2_vtx_IP ;

    TTree * _one_gamma_tree ;
    float _gamma_E;
    float _gamma_RL;
    float _gamma_vtx_IP ;


    TTree * _tree ;
    float _mu_startx ;
    float _mu_starty ;
    float _mu_startz ;
    float _mu_endx ;
    float _mu_endy ;
    float _mu_endz ;
    float _mu_mom;
    float _mu_len;
    float _mu_angle;
    float _mu_phi;
    float _mult;

    int _nshrs ;

    float _reco_E ;
    float _reco_start3D ;
    float _reco_startx ;
    float _reco_starty ;
    float _reco_startz ;
    float _reco_dot ;
    float _reco_dirx ;
    float _reco_diry ;
    float _reco_dirz ;

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
