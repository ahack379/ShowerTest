/**
 * \file GetInteractionInfo.h
 *
 * \ingroup CalcEfficiency
 * 
 * \brief Class def header for a class GetInteractionInfo
 *
 * @author ah673
 */

/** \addtogroup CalcEfficiency

    @{*/

#ifndef LARLITE_GETINTERACTIONINFO_H
#define LARLITE_GETINTERACTIONINFO_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class GetInteractionInfo
     User custom analysis class made by SHELL_USER_NAME
   */
  class GetInteractionInfo : public ana_base{
  
  public:

    /// Default constructor
    GetInteractionInfo(){ _name="GetInteractionInfo"; _fout=0;}

    /// Default destructor
    virtual ~GetInteractionInfo(){}

    /** IMPLEMENT in GetInteractionInfo.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GetInteractionInfo.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GetInteractionInfo.cc! 
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
