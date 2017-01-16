/**
 * \file EnergyFilter.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class EnergyFilter
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_ENERGYFILTER_H
#define LARLITE_ENERGYFILTER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class EnergyFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class EnergyFilter : public ana_base{
  
  public:

    /// Default constructor
    EnergyFilter(){ _name="EnergyFilter"; _fout=0;}

    /// Default destructor
    virtual ~EnergyFilter(){}

    /** IMPLEMENT in EnergyFilter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in EnergyFilter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in EnergyFilter.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:
    
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
