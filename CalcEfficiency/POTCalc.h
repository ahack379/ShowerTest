/**
 * \file POTCalc.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class POTCalc
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_POTCALC_H
#define LARLITE_POTCALC_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class POTCalc
     User custom analysis class made by SHELL_USER_NAME
   */
  class POTCalc : public ana_base{
  
  public:

    /// Default constructor
    POTCalc(){ _name="POTCalc"; _fout=0;}

    /// Default destructor
    virtual ~POTCalc(){}

    /** IMPLEMENT in POTCalc.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in POTCalc.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in POTCalc.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  int _event ;
  
  std::vector<int> _event_list ;
  std::multimap<float,float> _map ;

  float _tot_pot ;
    
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
