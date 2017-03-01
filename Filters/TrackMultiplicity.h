/**
 * \file TrackMultiplicity.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class TrackMultiplicity
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_TRACKMULTIPLICITY_H
#define LARLITE_TRACKMULTIPLICITY_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class TrackMultiplicity
     User custom analysis class made by SHELL_USER_NAME
   */
  class TrackMultiplicity : public ana_base{
  
  public:

    /// Default constructor
    TrackMultiplicity(){ _name="TrackMultiplicity"; _fout=0;}

    /// Default destructor
    virtual ~TrackMultiplicity(){}

    /** IMPLEMENT in TrackMultiplicity.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in TrackMultiplicity.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in TrackMultiplicity.cc! 
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
