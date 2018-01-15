/**
 * \file MisAlign.h
 *
 * \ingroup Misc
 * 
 * \brief Class def header for a class MisAlign
 *
 * @author ah673
 */

/** \addtogroup Misc

    @{*/

#ifndef LARLITE_MISALIGN_H
#define LARLITE_MISALIGN_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class MisAlign
     User custom analysis class made by SHELL_USER_NAME
   */
  class MisAlign : public ana_base{
  
  public:

    /// Default constructor
    MisAlign(){ _name="MisAlign"; _fout=0;}

    /// Default destructor
    virtual ~MisAlign(){}

    /** IMPLEMENT in MisAlign.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MisAlign.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MisAlign.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    int _event ;
    
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
