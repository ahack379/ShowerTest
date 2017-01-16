/**
 * \file HitThreshStudy.h
 %*
 * \ingroup HitThreshration
 * 
 * \brief Class def header for a class HitThreshStudy
 *
 * @author david caratelli
 */

/** \addtogroup HitThreshration

    @{*/

#ifndef LARLITE_HITTHRESHSTUDY_H
#define LARLITE_HITTHRESHSTUDY_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class HitThreshStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class HitThreshStudy : public ana_base{
  
  public:

    /// Default constructor
    HitThreshStudy(){ _name="HitThreshStudy"; _fout=0;}

    /// Default destructor
    virtual ~HitThreshStudy(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /// set hit producer
    void setHitProducer(std::string s) { _hit_producer = s; }
    void setClusProducer(std::string s) { _clus_producer = s; }

  protected:

    /// hit producer
    std::string _hit_producer;
    std::string _clus_producer;

    /// TTree
    TTree* _tree;
    double _reco_area, _reco_ampl;
    double _q; // charge form simch
    int _tick;
    int _pl;
    int _trk_size;
    int _wire, _chan;
    int _hit_multiplicity;
    double _x_st;
    double _y_st;
    double _z_st;
    double _px_st;
    double _py_st;
    double _pz_st;
    double _e_st;

    double _x_dp;
    double _y_dp;
    double _z_dp;
    double _px_dp;
    double _py_dp;
    double _pz_dp;
    double _e_dp;

    int _event ;
    int _clus;
    
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
