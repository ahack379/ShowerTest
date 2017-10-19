/**
 * \file NeutrinoMode.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class NeutrinoMode
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_NEUTRINOMODE_H
#define LARLITE_NEUTRINOMODE_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class NeutrinoMode
     User custom analysis class made by SHELL_USER_NAME
   */
  class NeutrinoMode : public ana_base{
  
  public:

    /// Default constructor
    NeutrinoMode(){ _name="NeutrinoMode"; _fout=0; _tree=0; }

    /// Default destructor
    virtual ~NeutrinoMode(){}

    /** IMPLEMENT in NeutrinoMode.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in NeutrinoMode.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in NeutrinoMode.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void clear() ;

  protected:

  std::vector<int> _event_list ;

  TTree * _tree ;
  int _event ;
  int _bkgd_id ;
  int _nu_mode;
    
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
