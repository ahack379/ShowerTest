/**
 * \file CCEff.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class CCEff
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_CCEFF_H
#define LARLITE_CCEFF_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class CCEff
     User custom analysis class made by SHELL_USER_NAME
   */
  class CCEff : public ana_base{
  
  public:

    /// Default constructor
    CCEff(){ _name="CCEff"; _fout=0;}

    /// Default destructor
    virtual ~CCEff(){}

    /** IMPLEMENT in CCEff.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CCEff.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CCEff.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  int _event ;
  int _signal ;

  int _n_cc;

  std::vector<int> _event_list ;
  std::multimap<float,float> _map;

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
