/**
 * \file RLEnergyCutStudy.h
 *
 * \ingroup FindCutValue
 * 
 * \brief Class def header for a class RLEnergyCutStudy
 *
 * @author ah673
 */

/** \addtogroup FindCutValue

    @{*/

#ifndef LARLITE_RLENERGYCUTSTUDY_H
#define LARLITE_RLENERGYCUTSTUDY_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class RLEnergyCutStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class RLEnergyCutStudy : public ana_base{
  
  public:

    /// Default constructor
    RLEnergyCutStudy(){ _name="RLEnergyCutStudy"; _fout = 0; _shower_tree = nullptr; }

    /// Default destructor
    virtual ~RLEnergyCutStudy(){}

    /** IMPLEMENT in RLEnergyCutStudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in RLEnergyCutStudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in RLEnergyCutStudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void clear();

  protected:

    int _event ;

    TTree * _shower_tree ;
    float _shower_e ;
    float _shower_rl ;
    bool _signal;
    
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
