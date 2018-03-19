/**
 * \file GetEventWeightsFromFullSample.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class GetEventWeightsFromFullSample
 *
 * @author ariana hackenburg 
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_GETEVENTWEIGHTSFROMFULLSAMPLE_H
#define LARLITE_GETEVENTWEIGHTSFROMFULLSAMPLE_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"

namespace larlite {
  /**
     \class GetEventWeightsFromFullSample
     User custom analysis class made by SHELL_USER_NAME
   */
  class GetEventWeightsFromFullSample : public ana_base{
  
  public:

    /// Default constructor
    GetEventWeightsFromFullSample() ; 

    /// Default destructor
    virtual ~GetEventWeightsFromFullSample(){}

    /** IMPLEMENT in GetEventWeightsFromFullSample.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GetEventWeightsFromFullSample.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GetEventWeightsFromFullSample.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetWeightProducer(std::string prod ) { _eventweight_label = prod; }

    void PrintOutput(bool doit ) { _print_output = doit; }

  protected:

    float _n_sig;
    
    /////////// My extra variables
    int _funcs ;
    std::vector<std::string> _genie_label_v ;
    std::vector<std::string> _unisim_label_v ;

    int _events ;
    std::string _eventweight_label ;
     
    float _cv;
    bool _signal ;

    std::vector<std::vector<float>>  _t_weights_by_universe ;

    std::map<std::string,int> _label_map ;

    bool _print_output ;

    float _all_evts_nominal;
    std::vector<float> _all_evts_m1 ;
    std::vector<float> _all_evts_p1 ;

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
