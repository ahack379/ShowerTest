/**
 * \file CCpi0Eff.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class CCpi0Eff
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_CCPI0EFF_H
#define LARLITE_CCPI0EFF_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class CCpi0Eff
     User custom analysis class made by SHELL_USER_NAME
   */
  class CCpi0Eff : public ana_base{
  
  public:

    /// Default constructor
    CCpi0Eff(){ _name="CCpi0Eff"; _fout=0;}

    /// Default destructor
    virtual ~CCpi0Eff(){}

    /** IMPLEMENT in CCpi0Eff.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CCpi0Eff.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CCpi0Eff.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  int _event ;
  int _signal ;

  std::vector<int> _event_list ;
  std::vector<int> _diff_list;
  std::vector<int> _fail_list;
  std::vector<int> _mes_v; 
  std::vector<int> _pi0_list ;

  int _multpi0 ;
  int _mesons ;
  int _leptons ;
  int _gammas;
  int _tot_ccpi0;
  int _nc_pi0;
    
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