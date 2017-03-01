/**
 * \file TrackStitcher.h
 *
 * \ingroup CalcEfficiency
 * 
 * \brief Class def header for a class TrackStitcher
 *
 * @author ariana Hackenburg
 */

/** \addtogroup CalcEfficiency

    @{*/

#ifndef LARLITE_TRACKSTITCHER_H
#define LARLITE_TRACKSTITCHER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class TrackStitcher
     User custom analysis class made by SHELL_USER_NAME
   */
  class TrackStitcher : public ana_base{
  
  public:

    /// Default constructor
    TrackStitcher(){ _name="TrackStitcher"; _fout=0; _tree=nullptr;} 

    /// Default destructor
    virtual ~TrackStitcher(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:

  std::vector<int> _trk_list;
  std::vector<int> _good_list;
  std::vector<int> _bad_list;

  int _event ;

  TTree * _tree ;
  std::vector<float> _length_v;
  int _track_mult ;
  int _nue_event ;
  int _nue_mult ;
  int _ccpi0_event ;
  int _ccpi0_mult ;
  float _nu_energy ;
  float _length ;
  int _pid ;
  float _vtx_diff ; 

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
