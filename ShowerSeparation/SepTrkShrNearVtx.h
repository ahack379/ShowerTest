/**
 * \file SepTrkShrNearVtx.h
 *
 * \ingroup ShowerSeparation
 * 
 * \brief Class def header for a class SepTrkShrNearVtx
 *
 * @author ah673
 */

/** \addtogroup ShowerSeparation

    @{*/

#ifndef LARLITE_SEPTRKSHRNEARVTX_H
#define LARLITE_SEPTRKSHRNEARVTX_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class SepTrkShrNearVtx
     User custom analysis class made by SHELL_USER_NAME
   */
  class SepTrkShrNearVtx : public ana_base{
  
  public:

    /// Default constructor
    SepTrkShrNearVtx(){ _name="SepTrkShrNearVtx"; _fout=0; _lin_tree=0;}

    /// Default destructor
    virtual ~SepTrkShrNearVtx(){}

    /** IMPLEMENT in SepTrkShrNearVtx.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SepTrkShrNearVtx.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SepTrkShrNearVtx.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void Clear();

  protected:

  int _event ; 

  TTree * _lin_tree;
  float _lin ;
  float _tll ;
  int _nhits ;
  bool _is_shower ;
  float _length ;

  std::vector<int> _event_list ; 
    
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
