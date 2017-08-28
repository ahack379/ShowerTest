/**
 * \file JustifyPi0Cuts.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class JustifyPi0Cuts
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_JUSTIFYPI0CUTS_H
#define LARLITE_JUSTIFYPI0CUTS_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/GeometryHelper.h"


namespace larlite {
  /**
     \class JustifyPi0Cuts
     User custom analysis class made by SHELL_USER_NAME
   */
  class JustifyPi0Cuts : public ana_base{
  
  public:

    /// Default constructor
    JustifyPi0Cuts(){ _name="JustifyPi0Cuts"; _fout=0; _pi0_selection=0; _tree=0; }

    /// Default destructor
    virtual ~JustifyPi0Cuts(){}

    /** IMPLEMENT in JustifyPi0Cuts.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in JustifyPi0Cuts.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in JustifyPi0Cuts.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void clear();

  protected:

    TTree * _pi0_selection;
    int _event;
    float _pi0_mass;
    float _pi0_mom;
    float _pi0_oangle;
    float _pi0_low_shrE;
    float _pi0_high_shrE;
    float _pi0_low_radL;
    float _pi0_high_radL;
    float _pi0_IP;

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
    int _bkgd_id;
    int _nshrs ;

    int _n_other ;
    int _n_cosmic ;
    int _n_cc1pi0 ; // This is the signal
    int _n_cc0pi0 ;
    int _n_nc1pi0 ; 
    int _n_nc0pi0 ;



    ::geoalgo::GeoAlgo _geoAlgo ;
    std::vector<int> _event_list ;    

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
