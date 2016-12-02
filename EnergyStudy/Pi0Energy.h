/**
 * \file Pi0Energy.h
 *
 * \ingroup EnergyStudy
 * 
 * \brief Class def header for a class Pi0Energy
 *
 * @author ah673
 */

/** \addtogroup EnergyStudy

    @{*/

#ifndef LARLITE_PI0ENERGY_H
#define LARLITE_PI0ENERGY_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class Pi0Energy
     User custom analysis class made by SHELL_USER_NAME
   */
  class Pi0Energy : public ana_base{
  
  public:

    /// Default constructor
    Pi0Energy(){ _name="Pi0Energy"; _fout=0; _energy_tree=0;}

    /// Default destructor
    virtual ~Pi0Energy(){}

    /** IMPLEMENT in Pi0Energy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in Pi0Energy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in Pi0Energy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void Clear() ;

  protected:

  TTree * _energy_tree ;
  float _true_e ;
  float _reco_e ;
  int _n_true_pi0;

  float _event ;
  std::vector<int> _event_list ;

  geoalgo::GeoAlgo _geoAlgo ;
    
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
