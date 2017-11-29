/**
 * \file FluxXSecErrorsSelected.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class FluxXSecErrorsSelected
 *
 * @author ariana hackenburg 
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_FLUXXSECERRORSSELECTED_H
#define LARLITE_FLUXXSECERRORSSELECTED_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"

namespace larlite {
  /**
     \class FluxXSecErrorsSelected
     User custom analysis class made by SHELL_USER_NAME
   */
  class FluxXSecErrorsSelected : public ana_base{
  
  public:

    /// Default constructor
    FluxXSecErrorsSelected() ; 

    /// Default destructor
    virtual ~FluxXSecErrorsSelected(){}

    /** IMPLEMENT in FluxXSecErrorsSelected.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FluxXSecErrorsSelected.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FluxXSecErrorsSelected.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetWeightProducer(std::string prod ) { _event_producer = prod; }

  protected:

    float _tot_pot ;
    const larutil::Geometry * fGeometry ;
    
    /////////// My extra variables
    float _xsec_mom_truth ;
    float _xsec_theta_truth ;
    std::vector<float> _weight_v ;

    float _sel_evts_nominal;
    std::vector<float> _sel_evts_m1 ;
    std::vector<float> _sel_evts_p1 ;
    std::vector<float> _sel_evts_nom ;

    float _bkgd_evts_nominal;
    std::vector<float> _bkgd_evts_m1 ;
    std::vector<float> _bkgd_evts_p1 ;

    std::vector<std::string> _genie_label_v ;

    int _events ;

    std::string _event_producer ;
     
    TTree * _tree ;
    float _cv;
    std::vector<float> _up ;
    std::vector<float> _down;
    bool _signal ;

    std::vector<std::vector<float>>  _s_weights_by_universe ;
    std::vector<std::vector<float>>  _b_weights_by_universe ;

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
