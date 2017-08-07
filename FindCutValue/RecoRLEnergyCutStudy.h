/**
 * \file RecoRLEnergyCutStudy.h
 *
 * \ingroup FindCutValue
 * 
 * \brief Class def header for a class RecoRLEnergyCutStudy
 *
 * @author ah673
 */

/** \addtogroup FindCutValue

    @{*/

#ifndef LARecoRLITE_RECORLENERGYCUTSTUDY_H
#define LARecoRLITE_RECORLENERGYCUTSTUDY_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class RecoRLEnergyCutStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class RecoRLEnergyCutStudy : public ana_base{
  
  public:

    /// Default constructor
    RecoRLEnergyCutStudy(){ _name="RecoRLEnergyCutStudy"; _fout = 0; _shower_tree = nullptr; _hit_tree = 0;}

    /// Default destructor
    virtual ~RecoRLEnergyCutStudy(){}

    /** IMPLEMENT in RecoRLEnergyCutStudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in RecoRLEnergyCutStudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in RecoRLEnergyCutStudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void clear();

  protected:


    TTree * _shower_tree ;
    int _event ;
    float _shower_e ;
    float _shower_rl ;
    float _shower_angle ;
    float _shower_oangle ;
    float _mu_angle ;
    bool _signal ;
    bool _nshrs;

    TTree * _hit_tree ;
    std::vector<float> _gaus_hits_per_r;
    std::vector<float> _shr_hits_per_r;
    std::vector<float> _radii; 
    std::vector<float> _hit_ratio ;

    ::geoalgo::GeoAlgo _geoAlgo ;

    
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
