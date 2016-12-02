/**
 * \file FullDistribs.h
 *
 * \ingroup EnergyStudy
 * 
 * \brief Class def header for a class FullDistribs
 *
 * @author ah673
 */

/** \addtogroup EnergyStudy

    @{*/

#ifndef LARLITE_FULLDISTRIBS_H
#define LARLITE_FULLDISTRIBS_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class FullDistribs
     User custom analysis class made by SHELL_USER_NAME
   */
  class FullDistribs : public ana_base{
  
  public:

    /// Default constructor
    FullDistribs(){ _name="FullDistribs"; _fout=0; _energy_tree=0; _gamma_tree=0;}

    /// Default destructor
    virtual ~FullDistribs(){}

    /** IMPLEMENT in FullDistribs.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FullDistribs.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FullDistribs.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void Clear() ;

  protected:

  TTree * _energy_tree ;
  float _true_pi0_e ;
  float _true_angle ;
  float _true_asym ;
  std::vector<float> _true_vtx ;
  std::vector<float> _reco_vtx ;

  int _n_true_pi0;

  float _event ;

  TTree * _gamma_tree ;
  float _gamma_e;
  float _rad_l ;

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
