/**
 * \file SingleShowerPi0Cuts.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class SingleShowerPi0Cuts
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_SINGLESHOWERPI0CUTS_H
#define LARLITE_SINGLESHOWERPI0CUTS_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/GeometryHelper.h"


namespace larlite {
  /**
     \class SingleShowerPi0Cuts
     User custom analysis class made by SHELL_USER_NAME
   */
  class SingleShowerPi0Cuts : public ana_base{
  
  public:

    /// Default constructor
    SingleShowerPi0Cuts(){ _name="SingleShowerPi0Cuts"; _fout=0; _pi0_selection=0; }

    /// Default destructor
    virtual ~SingleShowerPi0Cuts(){}

    /** IMPLEMENT in SingleShowerPi0Cuts.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SingleShowerPi0Cuts.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SingleShowerPi0Cuts.cc! 
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
    float _mu_mom;
    float _mu_angle;

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
