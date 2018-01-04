/**
 * \file FluxXSecMultWeights.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class FluxXSecMultWeights
 *
 * @author ariana hackenburg 
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_FLUXXSECMULTWEIGHTS_H
#define LARLITE_FLUXXSECMULTWEIGHTS_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"

namespace larlite {
  /**
     \class FluxXSecMultWeights
     User custom analysis class made by SHELL_USER_NAME
   */
  class FluxXSecMultWeights : public ana_base{
  
  public:

    /// Default constructor
    FluxXSecMultWeights() ; 

    /// Default destructor
    virtual ~FluxXSecMultWeights(){}

    /** IMPLEMENT in FluxXSecMultWeights.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FluxXSecMultWeights.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FluxXSecMultWeights.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetWeightProducer(std::string prod ) { _event_producer = prod; }

    void SetNominal ( float nom ) { _N = nom ; }

    void SetNominalXsec ( float nom ) { _N_xsec = nom ; }

  protected:

    float _tot_pot ;
    const larutil::Geometry * fGeometry ;
    
    /////////// My extra variables
    std::vector<float> _weight_v ;
    float _sel_evts_nominal;
    float _bkgd_evts_nominal;
    std::vector<std::string> _genie_label_v ;
    std::vector<std::string> _unisim_label_v ;

    int _events ;
    std::string _event_producer ;
     
    TTree * _tree ;
    bool _signal ;

    std::vector<std::vector<float>>  _s_weights_by_universe ;
    std::vector<std::vector<float>>  _b_weights_by_universe ;
    std::vector<std::vector<float>>  _t_weights_by_universe ;
    std::vector<std::vector<float>>  _flux_by_universe ; 

    TTree * _univ;
    std::vector<std::vector<float>> _xsec_v ;
    std::vector<std::vector<float>> _perc_v;

    std::map<std::string,int> _label_map ;


    float _N ;
    float _N_xsec ;

    int _bkgd_id ;
    std::map<int,int> _mc_hit_map ;
    float _mu_purity ;
    float _mu_complete;

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
