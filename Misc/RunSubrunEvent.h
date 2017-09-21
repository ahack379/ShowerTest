/**
 * \file RunSubrunEvent.h
 *
 * \ingroup Misc
 * 
 * \brief Class def header for a class RunSubrunEvent
 *
 * @author ah673
 */

/** \addtogroup Misc

    @{*/

#ifndef LARLITE_RUNSUBRUNEVENT_H
#define LARLITE_RUNSUBRUNEVENT_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class RunSubrunEvent
     User custom analysis class made by SHELL_USER_NAME
   */
  class RunSubrunEvent : public ana_base{
  
  public:

    /// Default constructor
    RunSubrunEvent(){ _name="RunSubrunEvent"; _fout=0;}

    /// Default destructor
    virtual ~RunSubrunEvent(){}

    /** IMPLEMENT in RunSubrunEvent.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in RunSubrunEvent.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in RunSubrunEvent.cc! 
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
