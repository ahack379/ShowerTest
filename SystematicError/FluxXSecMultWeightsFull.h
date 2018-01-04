/**
 * \file FluxXSecMultWeightsFull.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class FluxXSecMultWeightsFull
 *
 * @author ariana hackenburg 
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_FLUXXSECMULTWEIGHTSFULL_H
#define LARLITE_FLUXXSECMULTWEIGHTSFULL_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"

namespace larlite {
  /**
     \class FluxXSecMultWeightsFull
     User custom analysis class made by SHELL_USER_NAME
   */
  class FluxXSecMultWeightsFull : public ana_base{
  
  public:

    /// Default constructor
    FluxXSecMultWeightsFull() ; 

    /// Default destructor
    virtual ~FluxXSecMultWeightsFull(){}

    /** IMPLEMENT in FluxXSecMultWeightsFull.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FluxXSecMultWeightsFull.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FluxXSecMultWeightsFull.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetWeightProducer(std::string prod ) { _event_producer = prod; }

    void SetNominal ( float nom ) { _N = nom ; }

    void SetNominalXsec ( float nom ) { _N_xsec = nom ; }

  protected:

    float _n_sig;
    const larutil::Geometry * fGeometry ;
    
    /////////// My extra variables
    std::vector<std::string> _genie_label_v ;
    std::vector<std::string> _unisim_label_v ;

    int _events ;
    std::string _event_producer ;
     
    float _cv;
    bool _signal ;

    std::vector<std::vector<float>>  _t_weights_by_universe ;

    std::map<std::string,int> _label_map ;

    float _N ;
    float _N_xsec ;

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
