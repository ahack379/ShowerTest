/**
 * \file TrackMult.h
 *
 * \ingroup CalcEfficiency
 * 
 * \brief Class def header for a class TrackMult
 *
 * @author ariana Hackenburg
 */

/** \addtogroup CalcEfficiency

    @{*/

#ifndef LARLITE_TRACKMULT_H
#define LARLITE_TRACKMULT_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class TrackMult
     User custom analysis class made by SHELL_USER_NAME
   */
  class TrackMult : public ana_base{
  
  public:

    /// Default constructor
    TrackMult(){ _name="TrackMult"; _fout=0; _tree=nullptr; _ratio_cut = 0.; _radius = 40. ;}

    /// Default destructor
    virtual ~TrackMult(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetRatioCut( float r_cut ) { _ratio_cut = r_cut ; }

    void SetRadius( float rad ) { _radius = rad ; }

    bool IsSignal(int interactionType) ;

  protected:

  float _radius;          // Radius around vertex in which we'll check for hits
  float _ratio_cut ;      // Ratio cut on Shower-like : Gaus hits within radius
  int _good_tags ;        // Number of good events passing the filter
  int _bad_tags ;         // Number of bad events passing the filter
  int _total_sig_events;  // Number of total sig events for eff calculation

  std::vector<int> _trk_list;
  std::vector<int> _good_list;
  std::vector<int> _bad_list;

  int _event ;

  TTree * _tree ;
  std::vector<float> _length_v;
  int _track_mult ;
  bool _nue_event ;
  int _nue_mult ;
  bool _ccpi0_event ;
  int _ccpi0_mult ;
  float _nu_energy ;
  float _length ;
  int _pid ;



    
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
