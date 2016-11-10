/**
 * \file SeparateCCNC.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class SeparateCCNC
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_SEPARATECCNC_H
#define LARLITE_SEPARATECCNC_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class SeparateCCNC
     User custom analysis class made by SHELL_USER_NAME
   */
  class SeparateCCNC : public ana_base{
  
  public:

    /// Default constructor
    SeparateCCNC(){ _name="SeparateCCNC"; _fout=0; _getNC = false ;}

    /// Default destructor
    virtual ~SeparateCCNC(){}

    /** IMPLEMENT in SeparateCCNC.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SeparateCCNC.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SeparateCCNC.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void GetNC(bool get){ _getNC = get ; }

  protected:

  bool _getNC;
  int _signal ;

  int _events ;
    
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
