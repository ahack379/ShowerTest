/**
 * \file SignalEff.h
 *
 * \ingroup CalcEfficiency
 * 
 * \brief Class def header for a class SignalEff
 *
 * @author ah673
 */

/** \addtogroup CalcEfficiency

    @{*/

#ifndef LARLITE_SIGNALEFF_H
#define LARLITE_SIGNALEFF_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class SignalEff
     User custom analysis class made by SHELL_USER_NAME
   */
  class SignalEff : public ana_base{
  
  public:

    /// Default constructor
    SignalEff(){ _name="SignalEff"; _fout=0; _expected=0; _n_cosbkgd_evts=0; }

    /// Default destructor
    virtual ~SignalEff(){}

    /** IMPLEMENT in SignalEff.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SignalEff.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SignalEff.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    /// For now this will be a parameter passed from run time script
    void SetEventList(std::vector<int> sel_ev_list){ _sel_ev_list = sel_ev_list; }

    bool IsSelected(std::vector<int>, const int& event);

    void SetExpectedEvts(int expected){ _expected = expected; }

    void SetNCosmicBkgds(bool n_bkds ){ _n_cosbkgd_evts = n_bkds; }

  protected:

  std::vector<int> _sel_ev_list ;
  int _event ;

  int _sel_good ;
  int _sel_misid ;
  int _sel_tot ;  
  int _expected ;

  bool _n_cosbkgd_evts ;
    
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
