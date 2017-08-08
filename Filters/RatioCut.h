/**
 * \file RatioCut.h
 *
 * \ingroup SumCharge
 * 
 * \brief Class def header for a class RatioCut
 *
 * @author ah673
 */

/** \addtogroup SumCharge

    @{*/

#ifndef LARLITE_RATIOCUT_H
#define LARLITE_RATIOCUT_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class RatioCut
     User custom analysis class made by SHELL_USER_NAME
   */
  class RatioCut : public ana_base{
  
  public:

    /// Default constructor
    RatioCut(){ _name="RatioCut"; _fout=0; _tree=nullptr; _hits_tot=0.; _hits_in_rad=0; 
                  _ratio_cut=0.2; _radius = 60.; _vtx_producer="numuCC_vertex"; _chain_modules=false; }

    /// Default destructor
    virtual ~RatioCut(){}

    /** IMPLEMENT in RatioCut.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in RatioCut.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in RatioCut.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetRatioCut(float ratio ) { _ratio_cut = ratio ; } 

    void SetRadius(float rad ) { _radius = rad ; } 

    void SetVtxProducer(std::string vtx){ _vtx_producer = vtx ; } 

    void UseChainedModules(bool doit ) { _chain_modules = doit; }

  protected:

  float _ratio_cut ;
  float _radius ;
  std::string _vtx_producer ;
  bool _chain_modules ;

  TTree * _tree ;
  float _hits_tot;
  float _hits_in_rad; 
  float _hits_in_rad_g;

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
