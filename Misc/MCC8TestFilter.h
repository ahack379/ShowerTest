/**
 * \file MCC8TestFilter.h
 *
 * \ingroup Misc
 * 
 * \brief Class def header for a class MCC8TestFilter
 *
 * @author ah673
 */

/** \addtogroup Misc

    @{*/

#ifndef LARLITE_MCC8TESTFILTER_H
#define LARLITE_MCC8TESTFILTER_H

#include "Analysis/ana_base.h"
#include <map>

namespace larlite {
  /**
     \class MCC8TestFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class MCC8TestFilter : public ana_base{
  
  public:

    /// Default constructor
    MCC8TestFilter(){ _name="MCC8TestFilter"; _fout=0;}

    /// Default destructor
    virtual ~MCC8TestFilter(){}

    /** IMPLEMENT in MCC8TestFilter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MCC8TestFilter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MCC8TestFilter.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  std::map<int,std::pair<int,int>> _evt_m ;
    
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
