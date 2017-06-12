/**
 * \file ShowerCompleteness.h
 *
 * \ingroup Signal_Quality
 * 
 * \brief Class def header for a class ShowerCompleteness
 *
 * @author ah673
 */

/** \addtogroup Signal_Quality

    @{*/

#ifndef LARLITE_SHOWERCOMPLETENESS_H
#define LARLITE_SHOWERCOMPLETENESS_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class ShowerCompleteness
     User custom analysis class made by SHELL_USER_NAME
   */
  class ShowerCompleteness : public ana_base{
  
  public:

    /// Default constructor
    ShowerCompleteness(){ _name="ShowerCompleteness"; _fout=0; _tree = 0;}

    /// Default destructor
    virtual ~ShowerCompleteness(){}

    /** IMPLEMENT in ShowerCompleteness.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ShowerCompleteness.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ShowerCompleteness.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  // Multimap of plane -> hit index
  std::map<int,int> _mc_hit_map ;

  // First vector is plane, second vector is per cluster
  std::vector<std::vector<float>> _purity_v ;
  std::vector<std::vector<float>> _complete_v;

  TTree * _tree ;

  int _p0;
  int _p1;
  int _p2 ;

  int _tot_clus ;
  int _event ;
    
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
