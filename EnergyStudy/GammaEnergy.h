/**
 * \file GammaEnergy.h
 *
 * \ingroup EnergyStudy
 * 
 * \brief Class def header for a class GammaEnergy
 *
 * @author ah673
 */

/** \addtogroup EnergyStudy

    @{*/

#ifndef LARLITE_GAMMAENERGY_H
#define LARLITE_GAMMAENERGY_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class GammaEnergy
     User custom analysis class made by SHELL_USER_NAME
   */
  class GammaEnergy : public ana_base{
  
  public:

    /// Default constructor
    GammaEnergy(){ _name="GammaEnergy"; _fout=0; _tree=0;}

    /// Default destructor
    virtual ~GammaEnergy(){}

    /** IMPLEMENT in GammaEnergy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in GammaEnergy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in GammaEnergy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  TTree * _tree ; 
  int _event ;
  float _reco_e ;
  float _hit_reco_e ;
  float _mc_e ; 
  float _sum ;
  float _sum_adj ;
  float _sum_reco_int ;

    
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
