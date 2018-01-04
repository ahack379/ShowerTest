/**
 * \file DumbFilter.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class DumbFilter
 *
 * @author ahack379
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_DUMBFILTER_H
#define LARLITE_DUMBFILTER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class DumbFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class DumbFilter : public ana_base{
  
  public:

    /// Default constructor
    DumbFilter(){ _name="DumbFilter"; _fout=0;}

    /// Default destructor
    virtual ~DumbFilter(){}

    /** IMPLEMENT in DumbFilter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in DumbFilter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in DumbFilter.cc! 
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
