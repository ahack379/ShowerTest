/**
 * \file SysFindRoughEvents.h
 *
 * \ingroup Misc
 * 
 * \brief Class def header for a class SysFindRoughEvents
 *
 * @author ah673
 */

/** \addtogroup Misc

    @{*/

#ifndef LARLITE_SYSFINDROUGHEVENTS_H
#define LARLITE_SYSFINDROUGHEVENTS_H

#include <iostream>
#include <fstream>
#include <sstream>

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class SysFindRoughEvents
     User custom analysis class made by SHELL_USER_NAME
   */
  class SysFindRoughEvents : public ana_base{
  
  public:

    /// Default constructor
    SysFindRoughEvents(){ _name="SysFindRoughEvents"; _fout=0;}

    /// Default destructor
    virtual ~SysFindRoughEvents(){}

    /** IMPLEMENT in SysFindRoughEvents.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SysFindRoughEvents.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SysFindRoughEvents.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  std::map<float,float> _sel2_ext_m ;
  std::map<float,float> _pi0_ext_m ;

  std::ifstream _file;

    
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
