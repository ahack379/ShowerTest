/**
 * \file TwoShowerPi0Cuts.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class TwoShowerPi0Cuts
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_TWOSHOWERPI0CUTS_H
#define LARLITE_TWOSHOWERPI0CUTS_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/GeometryHelper.h"


namespace larlite {
  /**
     \class TwoShowerPi0Cuts
     User custom analysis class made by SHELL_USER_NAME
   */
  class TwoShowerPi0Cuts : public ana_base{
  
  public:

    /// Default constructor
    TwoShowerPi0Cuts(){ _name="TwoShowerPi0Cuts"; _fout=0; _pi0_selection=0; _chain_modules=false; }

    /// Default destructor
    virtual ~TwoShowerPi0Cuts(){}

    /** IMPLEMENT in TwoShowerPi0Cuts.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in TwoShowerPi0Cuts.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in TwoShowerPi0Cuts.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void clear();

    void UseChainedModules(bool doityouwont){ _chain_modules = doityouwont; } 

    void SetVtxProducer(std::string vtx_producer="numuCC_vertex" ){ _vtx_producer = vtx_producer; } 

  protected:

    bool _chain_modules ;

    TTree * _pi0_selection;
    int _event;
    float _pi0_mass;
    float _pi0_mom;
    float _pi0_oangle;
    float _pi0_low_shrE;
    float _pi0_high_shrE;
    float _pi0_low_radL;
    float _pi0_high_radL;

    float _mu_startx ;
    float _mu_starty ;
    float _mu_startz ;
    float _mu_endx ;
    float _mu_endy ;
    float _mu_endz ;
    float _mu_mom;
    float _mu_len;
    float _mu_angle;
    float _mu_phi;
    float _mult;

    ::geoalgo::GeoAlgo _geoAlgo ;
    std::vector<int> _event_list ;    
    std::vector<int> _one_shower_list ;

    std::string _vtx_producer ;

    int _more_than_one ;

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
