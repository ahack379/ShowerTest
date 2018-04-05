/**
 * \file EvtWeightFilter.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class EvtWeightFilter
 *
 * @author ah673
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_EVTWEIGHTFILTER_H
#define LARLITE_EVTWEIGHTFILTER_H

#include "Analysis/ana_base.h"
#include <fstream>

namespace larlite {
  /**
     \class EvtWeightFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class EvtWeightFilter : public ana_base{
  
  public:

    /// Default constructor
    EvtWeightFilter(){ _name="EvtWeightFilter"; _fout=0;}

    /// Default destructor
    virtual ~EvtWeightFilter(){}

    /** IMPLEMENT in EvtWeightFilter.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in EvtWeightFilter.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in EvtWeightFilter.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  std::ifstream _file;
  std::vector<std::multimap<float,float>> _pi0_map;
  std::vector<std::multimap<float,float>> _map_v;
  //std::multimap<float,float> _pi0_map;
  //
  int _event;
    
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
