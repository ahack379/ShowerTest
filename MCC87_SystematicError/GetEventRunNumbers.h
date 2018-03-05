/**
 * \file GetEventRunNumbers.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class GetEventRunNumbers
 *
 * @author ah673
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_GETEVENTRUNNUMBERS_H
#define LARLITE_GETEVENTRUNNUMBERS_H

#include "Analysis/ana_base.h"
#include "DataFormat/user_info.h"

namespace larlite {
  /**
     \class GetEventRunNumbers
     User custom analysis class made by SHELL_USER_NAME
   */
  class GetEventRunNumbers : public ana_base{
  
  public:

    /// Default constructor
    GetEventRunNumbers(){ _name="GetEventRunNumbers"; _fout=0; _ev_user=0; }

    /// Default destructor
    virtual ~GetEventRunNumbers(){}

    /** IMPLEMENT in GetEventRunNumbers.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GetEventRunNumbers.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GetEventRunNumbers.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  std::multimap<float,float> _map;

  int evt ;

  ::larlite::event_user* _ev_user ;
    
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
