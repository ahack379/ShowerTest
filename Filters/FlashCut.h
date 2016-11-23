/**
 * \file FlashCut.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class FlashCut
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_FLASHCUT_H
#define LARLITE_FLASHCUT_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class FlashCut
     User custom analysis class made by SHELL_USER_NAME
   */
  class FlashCut : public ana_base{
  
  public:

    /// Default constructor
    FlashCut(){ _name="FlashCut"; _fout=0;}

    /// Default destructor
    virtual ~FlashCut(){}

    /** IMPLEMENT in FlashCut.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FlashCut.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FlashCut.cc! 
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
