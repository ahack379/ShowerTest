/**
 * \file CrossSection.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class CrossSection
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_CROSSSECTION_H
#define LARLITE_CROSSSECTION_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class CrossSection
     User custom analysis class made by SHELL_USER_NAME
   */
  class CrossSection : public ana_base{
  
  public:

    /// Default constructor
    CrossSection(){ _name="CrossSection"; _fout=0;}

    /// Default destructor
    virtual ~CrossSection(){}

    /** IMPLEMENT in CrossSection.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CrossSection.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CrossSection.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  int _event ;
  int _signal;
  int _tot_event_in_AV ;
  int _n_numu ;
  int _n_nu_all ;

  double _tot_pot ;

  float _xmin ;
  float _xmax ;
  float _ymin ; 
  float _ymax ; 
  float _zmin ; 
  float _zmax ; 

  float _mean_e ;

  std::vector<int> _event_list ;
    
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
