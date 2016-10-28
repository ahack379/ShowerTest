/**
 * \file VtxDensity.h
 *
 * \ingroup SumCharge
 * 
 * \brief Class def header for a class VtxDensity
 *
 * @author ah673
 */

/** \addtogroup SumCharge

    @{*/

#ifndef LARLITE_VTXDENSITY_H
#define LARLITE_VTXDENSITY_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class VtxDensity
     User custom analysis class made by SHELL_USER_NAME
   */
  class VtxDensity : public ana_base{
  
  public:

    /// Default constructor
    VtxDensity(){ _name="VtxDensity"; _fout=0; _tree=nullptr; _hits_tot=0.; _hits_in_rad=0; 
                  _use_mcbnb_info = false; }

    /// Default destructor
    virtual ~VtxDensity(){}

    /** IMPLEMENT in VtxDensity.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in VtxDensity.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in VtxDensity.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void UseMCBNBInfo(bool setMC ) { _use_mcbnb_info = setMC ; }

  protected:

  TTree * _tree ;
  float _hits_tot;
  float _hits_in_rad; 
  float _hits_in_rad_g;
  bool _use_mcbnb_info ;

  int _event ;

  std::vector<float> _radii;
  std::vector<float> _density;
  std::vector<float> _hits_per_r;
    
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
