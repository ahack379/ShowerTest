/**
 * \file FluxXSecErrorsFull.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class FluxXSecErrorsFull
 *
 * @author ariana hackenburg 
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_FLUXXSECERRORSFULL_H
#define LARLITE_FLUXXSECERRORSFULL_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/Geometry.h"
#include <fstream>

namespace larlite {
  /**
     \class FluxXSecErrorsFull
     User custom analysis class made by SHELL_USER_NAME
   */
  class FluxXSecErrorsFull : public ana_base{
  
  public:

    /// Default constructor
    FluxXSecErrorsFull(); //{ _name="FluxXSecErrorsFull"; _fout=0;}

    /// Default destructor
    virtual ~FluxXSecErrorsFull(){}

    /** IMPLEMENT in FluxXSecErrorsFull.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in FluxXSecErrorsFull.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in FluxXSecErrorsFull.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    std::ifstream _file; //("list.txt",std::ifstream::in);

    float _tot_pot ;
    const larutil::Geometry * fGeometry ;

    std::multimap<float,float> _map ;
    
    /////////// My extra variables
    TTree * _tree ;
    float _xsec_mom_truth ;
    float _xsec_theta_truth ;
    std::vector<float> _weight_v ;

    TTree * _final_tree ;
    float _all_evts_nominal;
    std::vector<float> _all_evts_m1 ;
    std::vector<float> _all_evts_p1 ;
    std::vector<std::string> _genie_label_v ;

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
