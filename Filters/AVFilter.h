/**
 * \file AVFilter.h
 *
 * \ingroup CalcEfficiency
 * 
 * \brief Class def header for a class AVFilter
 *
 * @author ah673
 */

/** \addtogroup CalcEfficiency

    @{*/

#ifndef LARLITE_AVFILTER_H
#define LARLITE_AVFILTER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class AVFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class AVFilter : public ana_base{
  
  public:

    /// Default constructor
    AVFilter(){ _name="AVFilter"; _fout=0;}

    /// Default destructor
    virtual ~AVFilter(){}

    /** IMPLEMENT in AVFilter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in AVFilter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in AVFilter.cc! 
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
