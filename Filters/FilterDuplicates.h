/**
 * \file FilterDuplicates.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class FilterDuplicates
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_FILTERDUPLICATES_H
#define LARLITE_FILTERDUPLICATES_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class FilterDuplicates
     User custom analysis class made by SHELL_USER_NAME
   */
  class FilterDuplicates : public ana_base{
  
  public:

    /// Default constructor
    FilterDuplicates(){ _name="FilterDuplicates"; _fout=0;}

    /// Default destructor
    virtual ~FilterDuplicates(){}

    /** IMPLEMENT in FilterDuplicates.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FilterDuplicates.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FilterDuplicates.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  std::vector<std::multimap<float,float>> _map_v;

  int _event ;
  int _event_no_dup ;
    
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
