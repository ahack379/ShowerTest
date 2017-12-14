/**
 * \file TestSlimmedMCPart.h
 *
 * \ingroup Misc
 * 
 * \brief Class def header for a class TestSlimmedMCPart
 *
 * @author ah673
 */

/** \addtogroup Misc

    @{*/

#ifndef LARLITE_TESTSLIMMEDMCPART_H
#define LARLITE_TESTSLIMMEDMCPART_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class TestSlimmedMCPart
     User custom analysis class made by SHELL_USER_NAME
   */
  class TestSlimmedMCPart : public ana_base{
  
  public:

    /// Default constructor
    TestSlimmedMCPart(){ _name="TestSlimmedMCPart"; _fout=0;}

    /// Default destructor
    virtual ~TestSlimmedMCPart(){}

    /** IMPLEMENT in TestSlimmedMCPart.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in TestSlimmedMCPart.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in TestSlimmedMCPart.cc! 
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
