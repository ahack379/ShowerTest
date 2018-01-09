/**
 * \file ExtraDumbFilter.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class ExtraDumbFilter
 *
 * @author ahack379
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_EXTRADUMBFILTER_H
#define LARLITE_EXTRADUMBFILTER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class ExtraDumbFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class ExtraDumbFilter : public ana_base{
  
  public:

    /// Default constructor
    ExtraDumbFilter(){ _name="ExtraDumbFilter"; _fout=0;}

    /// Default destructor
    virtual ~ExtraDumbFilter(){}

    /** IMPLEMENT in ExtraDumbFilter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ExtraDumbFilter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ExtraDumbFilter.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetEvent( int event ) { _event = event ; }
    void SetSubrun( int subrun ) { _subrun = subrun ; }

  protected:
    
  int _event ;
  int _subrun ;

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
