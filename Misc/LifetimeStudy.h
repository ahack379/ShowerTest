/**
 * \file LifetimeStudy.h
 *
 * \ingroup Misc
 * 
 * \brief Class def header for a class LifetimeStudy
 *
 * @author ah673
 */

/** \addtogroup Misc

    @{*/

#ifndef LARLITE_LIFETIMESTUDY_H
#define LARLITE_LIFETIMESTUDY_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class LifetimeStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class LifetimeStudy : public ana_base{
  
  public:

    /// Default constructor
    LifetimeStudy(){ _name="LifetimeStudy"; _fout=0; _trk_tree = 0; _shr_tree = 0; _all_tree = 0; }

    /// Default destructor
    virtual ~LifetimeStudy(){}

    /** IMPLEMENT in LifetimeStudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in LifetimeStudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in LifetimeStudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TTree *_trk_tree;
    float _flash_t;
    float _flash_pe;
    float _hit_t;
    float _hit_charge;

    TTree *_shr_tree;
    float _hit_shr_t;
    float _hit_shr_charge;

    TTree * _all_tree;
    int _is_shower ;
    int _pdg;
    int _mom_pdg;
    
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
