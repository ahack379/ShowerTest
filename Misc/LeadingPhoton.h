/**
 * \file LeadingPhoton.h
 *
 * \ingroup Misc 
 * 
 * \brief Class def header for a class LeadingPhoton
 *
 * @author ah673
 */

/** \addtogroup Misc 

    @{*/

#ifndef LARLITE_LEADINGPHOTON_H
#define LARLITE_LEADINGPHOTON_H

#include "Analysis/ana_base.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class LeadingPhoton
     User custom analysis class made by SHELL_USER_NAME
   */
  class LeadingPhoton : public ana_base{
  
  public:

    /// Default constructor
    LeadingPhoton(){ _name="LeadingPhoton"; _fout=0; _tree=0; }

    /// Default destructor
    virtual ~LeadingPhoton(){}

    /** IMPLEMENT in LeadingPhoton.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in LeadingPhoton.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in LeadingPhoton.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  std::vector<int> _event_list ;
  
  int _one_shower_events ;
  int _event;

  TTree * _tree ;
  float _dot_prod ;
  float _reco_e ;
  float _mc_e ;

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
