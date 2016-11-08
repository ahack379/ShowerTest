/**
 * \file SeparateBNB.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class SeparateBNB
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_SEPARATEBNB_H
#define LARLITE_SEPARATEBNB_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class SeparateBNB
     User custom analysis class made by SHELL_USER_NAME
   */
  class SeparateBNB : public ana_base{
  
  public:

    /// Default constructor
    SeparateBNB(){ _name="SeparateBNB"; _fout=0; _get_shower_events = false ;}

    /// Default destructor
    virtual ~SeparateBNB(){}

    /** IMPLEMENT in SeparateBNB.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SeparateBNB.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SeparateBNB.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void GetShowerEvts(bool get){ _get_shower_events = get ; }

  protected:

  bool _get_shower_events; 
    
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
