/**
 * \file MCVariableStudy.h
 *
 * \ingroup EnergyStudy
 * 
 * \brief Class def header for a class MCVariableStudy
 *
 * @author ah673
 */

/** \addtogroup EnergyStudy

    @{*/

#ifndef LARLITE_MCVARIABLESTUDY_H
#define LARLITE_MCVARIABLESTUDY_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class MCVariableStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class MCVariableStudy : public ana_base{
  
  public:

    /// Default constructor
    MCVariableStudy(){ _name="MCVariableStudy"; _fout=0; _pi0_tree=0; _gamma_tree=0;}

    /// Default destructor
    virtual ~MCVariableStudy(){}

    /** IMPLEMENT in MCVariableStudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MCVariableStudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MCVariableStudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void Clear() ;

  protected:

  TTree * _pi0_tree ;
  int _run_id ;
  int _subrun_id ;
  int _event_id ;
  float _true_pi0_e ;
  float _true_2mcs_e ;
  float _true_angle ;
  float _true_asym ;
  float _true_gamma_e_max ;
  float _true_gamma_e_min ;
  float _true_RL_maxE ;
  float _true_RL_minE ;
  float _true_pi0_mom ;
  float _true_nu_e ;
  float _true_mu_mom ;
  float _true_mu_len ;
  float _true_mu_theta;
  float _true_mu_phi;

  float _reco_pi0_e ;
  int _n_true_pi0;

  float _event ;

  TTree * _gamma_tree ;
  float _true_gamma_e;
  float _reco_gamma_e;
  float _true_rad_l ;
  float _reco_rad_l ;
  float _true_reco_dot ;


  std::vector<int> _event_list ;

  geoalgo::GeoAlgo _geoAlgo ;

  int _n_signals ;
    
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
