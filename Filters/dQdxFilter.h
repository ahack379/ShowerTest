/**
 * \file dQdxFilter.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class dQdxFilter
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_DQDXFILTER_H
#define LARLITE_DQDXFILTER_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class dQdxFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class dQdxFilter : public ana_base{
  
  public:

    /// Default constructor
    dQdxFilter(){ _name="dQdxFilter"; _fout=0; _use_onbeam = false; _use_offbeam = false; _use_mcbnbcos = false;}

    /// Default destructor
    virtual ~dQdxFilter(){}

    /** IMPLEMENT in dQdxFilter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in dQdxFilter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in dQdxFilter.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void UseOn( bool doit=false ){_use_onbeam = doit; }

    void UseOff( bool doit=false ){_use_offbeam = doit; }

    void UseMC( bool doit=false ){_use_mcbnbcos = doit; }

  protected:

  std::vector<int> _onbeam_v;
  std::vector<int> _offbeam_v ;
  std::vector<int> _mcbnbcos_v ;

  std::vector<int> _ev_v;

  bool _use_onbeam ;
  bool _use_offbeam ;
  bool _use_mcbnbcos ;
    
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
