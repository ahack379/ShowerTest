/**
 * \file AngledQdxCuts.h
 *
 * \ingroup Filters
 * 
 * \brief Class def header for a class AngledQdxCuts
 *
 * @author ah673
 */

/** \addtogroup Filters

    @{*/

#ifndef LARLITE_ANGLEDQDXCUTS_H
#define LARLITE_ANGLEDQDXCUTS_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"

namespace larlite {
  /**
     \class AngledQdxCuts
     User custom analysis class made by SHELL_USER_NAME
   */
  class AngledQdxCuts : public ana_base{
  
  public:

    /// Default constructor
    AngledQdxCuts(){ _name="AngledQdxCuts"; _fout=0; _use_onbeam = false; _use_offbeam = false; _use_mcbnbcos = false;}

    /// Default destructor
    virtual ~AngledQdxCuts(){}

    /** IMPLEMENT in AngledQdxCuts.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in AngledQdxCuts.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in AngledQdxCuts.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    double Median(std::vector<double> in) ;
    double TrunMean(std::vector<double> in) ;
    double MaxDeflection(larlite::track t );

    void UseMC( bool doit ){ _use_mcbnbcos = doit; }

  protected:

  std::vector<int> _onbeam_v;
  std::vector<int> _offbeam_v ;
  std::vector<int> _mcbnbcos_v ;

  std::vector<int> _ev_v;

  bool _use_onbeam ;
  bool _use_offbeam ;
  bool _use_mcbnbcos ;

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
