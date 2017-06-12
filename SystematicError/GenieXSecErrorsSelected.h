/**
 * \file GenieXSecErrorsSelected.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class GenieXSecErrorsSelected
 *
 * @author ariana hackenburg 
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_GENIEXSECERRORSSELECTED_H
#define LARLITE_GENIEXSECERRORSSELECTED_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"

namespace larlite {
  /**
     \class GenieXSecErrorsSelected
     User custom analysis class made by SHELL_USER_NAME
   */
  class GenieXSecErrorsSelected : public ana_base{
  
  public:

    /// Default constructor
    GenieXSecErrorsSelected(); //{ _name="GenieXSecErrorsSelected"; _fout=0;}

    /// Default destructor
    virtual ~GenieXSecErrorsSelected(){}

    /** IMPLEMENT in GenieXSecErrorsSelected.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GenieXSecErrorsSelected.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GenieXSecErrorsSelected.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    float _tot_pot ;
    const larutil::Geometry * fGeometry ;
    
    /////////// My extra variables
    TTree * _tree ;
    float _xsec_mom_truth ;
    float _xsec_theta_truth ;
    std::vector<float> _weight_v ;

    TTree * _final_tree ;
    float _sel_evts_nominal;
    std::vector<float> _sel_evts_m1 ;
    std::vector<float> _sel_evts_p1 ;

    float _bkgd_evts_nominal;
    std::vector<float> _bkgd_evts_m1 ;
    std::vector<float> _bkgd_evts_p1 ;

    std::vector<std::string> _genie_label_v ;

    int _events ;

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
