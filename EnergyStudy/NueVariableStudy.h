/**
 * \file NueVariableStudy.h
 *
 * \ingroup EnergyStudy
 * 
 * \brief Class def header for a class NueVariableStudy
 *
 * @author ah673
 */

/** \addtogroup EnergyStudy

    @{*/

#ifndef LARLITE_NUEVARIABLESTUDY_H
#define LARLITE_NUEVARIABLESTUDY_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class NueVariableStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class NueVariableStudy : public ana_base{
  
  public:

    /// Default constructor
    NueVariableStudy(){ _name="NueVariableStudy"; _fout=0; _nue_selection=0; _pot_tree=0;}

    /// Default destructor
    virtual ~NueVariableStudy(){}

    /** IMPLEMENT in NueVariableStudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in NueVariableStudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in NueVariableStudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  TTree * _nue_selection ;
  float _energy; 
  float _mc_xvtx ;
  float _mc_yvtx ;
  float _mc_zvtx ;
  float _rc_xvtx ;
  float _rc_yvtx ;
  float _rc_zvtx ;
  float _vtx_dist;
  int _pass ;
  bool _nc ;

  float _event ;
  std::vector<int> _event_list ;

  TTree * _pot_tree ;
  float _pottotl ;
  float _evttotl ;
    
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
