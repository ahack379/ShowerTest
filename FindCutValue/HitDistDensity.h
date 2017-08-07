/**
 * \file HitDistDensity.h
 *
 * \ingroup SumCharge
 * 
 * \brief Class def header for a class HitDistDensity
 *
 * @author ah673
 */

/** \addtogroup SumCharge

    @{*/

#ifndef LARLITE_HITDISTDENSITY_H
#define LARLITE_HITDISTDENSITY_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class HitDistDensity
     User custom analysis class made by SHELL_USER_NAME
   */
  class HitDistDensity : public ana_base{
  
  public:

    /// Default constructor
    HitDistDensity(){ _name="HitDistDensity"; _fout=0; _tree=nullptr; _hits_tot=0; _hits_in_rad=0; _get_pi0s = false; }

    /// Default destructor
    virtual ~HitDistDensity(){}

    /** IMPLEMENT in HitDistDensity.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in HitDistDensity.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in HitDistDensity.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void GetPi0s( bool getem) { _get_pi0s = getem ; }

  protected:

  TTree * _tree ;
  int _hits_tot;
  float _hits_in_rad; 
  float _hits_in_rad_g;

  int plane_to_use = 2;

  int _event ;
  int _keep ;

  std::vector<float> _radii;
  std::vector<float> _density;
  std::vector<float> _hits_per_r;
  std::vector<float> _shr_hits_per_r;
  std::vector<float> _gaus_hits_per_r;

  float _vtx_x ;
  float _vtx_y ;
  float _vtx_z ;

  bool _get_pi0s;

  int _shr_hits_tot;
    
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
