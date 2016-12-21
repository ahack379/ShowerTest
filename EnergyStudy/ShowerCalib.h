/**
 * \file ShowerCalib.h
 *
 * \ingroup EnergyStudy
 * 
 * \brief Class def header for a class ShowerCalib
 *
 * @author ah673
 */

/** \addtogroup EnergyStudy

    @{*/

#ifndef LARLITE_SHOWERCALIB_H
#define LARLITE_SHOWERCALIB_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class ShowerCalib
     User custom analysis class made by SHELL_USER_NAME
   */
  class ShowerCalib : public ana_base{
  
  public:

    /// Default constructor
    ShowerCalib(){ _name="ShowerCalib"; _fout=0; _pi0_tree=0; _gamma_tree=0;}

    /// Default destructor
    virtual ~ShowerCalib(){}

    /** IMPLEMENT in ShowerCalib.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ShowerCalib.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ShowerCalib.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void Clear() ;

  protected:

  TTree * _pi0_tree ;
  float _true_pi0_e ;
  float _true_2mcs_e ;
  float _true_angle ;
  float _true_asym ;

  float _reco_pi0_e ;
  int _n_true_pi0;

  float _event ;

  TTree * _gamma_tree ;
  float _true_gamma_e;
  float _reco_gamma_e;
  float _true_adj_gamma_e;
  float _reco_adj_gamma_e;
  float _true_rad_l ;
  float _reco_rad_l ;
  float _true_reco_dot ;

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
