/**
 * \file JustifyPi0Cuts.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class JustifyPi0Cuts
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_JUSTIFYCUTS_H
#define LARLITE_JUSTIFYCUTS_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"


namespace larlite {
  /**
     \class JustifyPi0Cuts
     User custom analysis class made by SHELL_USER_NAME
   */
  class JustifyPi0Cuts : public ana_base{
  
  public:

    /// Default constructor
    JustifyPi0Cuts(){ _name="JustifyPi0Cuts"; _fout=0; _gamma_tree=0; _tree=0; _one_gamma_tree = 0; _compare_tree = 0; }

    /// Default destructor
    virtual ~JustifyPi0Cuts(){}

    /** IMPLEMENT in JustifyPi0Cuts.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in JustifyPi0Cuts.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in JustifyPi0Cuts.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void clear();

    void clearCompare();

  protected:

    TTree * _gamma_tree ;
    int _event;
    int _event_type; // signal 0, bkgd with pi0 1, bkgd with no pi0 2
    float _pi0_mass;
    float _pi0_mom;
    float _gamma_oangle;
    float _gamma_low_E;
    float _gamma_high_E;
    float _gamma_low_RL;
    float _gamma_high_RL;
    float _gamma_IP;
    bool _gamma_low_matched ;
    bool _gamma_high_matched ;
    float _gamma1_vtx_IP ;
    float _gamma2_vtx_IP ;

    float _gamma_startx ;
    float _gamma_starty ;
    float _gamma_startz ;


    int _nu_pdg ;
    bool _isCC ;
    bool _found_pi0 ; 
    int _n_nu_origin_pi0 ; 

    float _gamma_low_purity ;
    float _gamma_low_complete ;
    float _gamma_high_purity ;
    float _gamma_high_complete ;

    float _mc_low_dirx ;
    float _mc_low_diry ;
    float _mc_low_dirz ;
    float _mc_high_dirx ;
    float _mc_high_diry ;
    float _mc_high_dirz ;
    float _gamma_low_dirx ;
    float _gamma_low_diry ;
    float _gamma_low_dirz ;
    float _gamma_high_dirx ;
    float _gamma_high_diry ;
    float _gamma_high_dirz ;


    float _mc_low_startx ;
    float _mc_low_starty ;
    float _mc_low_startz ;
    float _mc_high_startx ;
    float _mc_high_starty ;
    float _mc_high_startz ;
    float _gamma_low_startx ;
    float _gamma_low_starty ;
    float _gamma_low_startz ;
    float _gamma_high_startx ;
    float _gamma_high_starty ;
    float _gamma_high_startz ;

    float _res_high ;
    float _res_low ;

    bool _pi0_origin ;
    bool _pi0_type;

    TTree * _one_gamma_tree ;
    float _gamma_E;
    float _gamma_RL;
    float _gamma_vtx_IP ;
    bool _gamma_matched ;

    float _gamma_purity ;
    float _gamma_complete;
    float _mc_dirx ;
    float _mc_diry ;
    float _mc_dirz ;

    float _mc_startx ;
    float _mc_starty ;
    float _mc_startz ;
    float _gamma_dirx ;
    float _gamma_diry ;
    float _gamma_dirz ;
    float _res ;

    int _gamma_origin ;
    int _gamma_type;


    TTree * _tree ;
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
    int _bkgd_id;
    int _nshrs ;

    TTree * _compare_tree ;
    float _reco_E ;
    float _reco_start3D ;
    float _reco_startx ;
    float _reco_starty ;
    float _reco_startz ;
    float _reco_dot ;
    float _reco_dirx ;
    float _reco_diry ;
    float _reco_dirz ;

    float _mc_E ;
    //float _mc_startx ;
    //float _mc_starty ;
    //float _mc_startz ;
    //float _mc_dirx ;
    //float _mc_diry ;
    //float _mc_dirz ;


    int _n_other ;
    int _n_cosmic ;
    int _n_cc1pi0 ; // This is the signal
    int _n_cc0pi0 ;
    int _n_nc1pi0 ; 
    int _n_nc0pi0 ;


    ::geoalgo::GeoAlgo _geoAlgo ;
    double _time2cm;
    larutil::SpaceChargeMicroBooNE *_SCE;

    std::map<int,int> _mc_hit_map ;



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
