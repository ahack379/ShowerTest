/**
 * \file FilterHRNoMuon.h
 *
 * \ingroup CalcEfficiency
 * 
 * \brief Class def header for a class FilterHRNoMuon
 *
 * @author ariana Hackenburg
 */

/** \addtogroup CalcEfficiency

    @{*/

#ifndef LARLITE_FILTERHRNOMUON_H
#define LARLITE_FILTERHRNOMUON_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class FilterHRNoMuon
     User custom analysis class made by SHELL_USER_NAME
   */
  class FilterHRNoMuon : public ana_base{
  
  public:

    /// Default constructor
    FilterHRNoMuon(){ _name="FilterHRNoMuon"; _fout=0; _tree=nullptr; _ratio_cut = 0.; _radius = 40. ;}

    /// Default destructor
    virtual ~FilterHRNoMuon(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetRatioCut( float r_cut ) { _ratio_cut = r_cut ; }

    void SetRadius( float rad ) { _radius = rad ; }

    bool IsSignal(int interactionType) ;

  protected:

  float _radius;          // Radius around vertex in which we'll check for hits
  float _hits_in_rad_g;   // Gaus hits in radius
  float _hits_in_rad;     // Shower-like hits in radius
  float _hits_in_rad_ass; // Association track hits in radius
  float _ratio_cut ;      // Ratio cut on Shower-like : Gaus hits within radius
  int _good_tags ;        // Number of good events passing the filter
  int _bad_tags ;         // Number of bad events passing the filter
  int _total_sig_events;  // Number of total sig events for eff calculation

  int _event ;
  std::vector<float> _hit_ratio_v ;
  float _hit_ratio ;
  TTree * _tree ;

    
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
