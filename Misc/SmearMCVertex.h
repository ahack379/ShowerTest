/**
 * \file SmearMCVertex.h
 *
 * \ingroup Misc
 * 
 * \brief Class def header for a class SmearMCVertex
 *
 * @author ah673
 */

/** \addtogroup Misc

    @{*/

#ifndef LARLITE_SMEARMCVERTEX_H
#define LARLITE_SMEARMCVERTEX_H

#include "Analysis/ana_base.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "TH1D.h"

namespace larlite {
  /**
     \class SmearMCVertex
     User custom analysis class made by SHELL_USER_NAME
   */
  class SmearMCVertex : public ana_base{
  
  public:

    /// Default constructor
    SmearMCVertex(){ _name="SmearMCVertex"; _fout=0; _hist_x=0; _hist_y=0; _hist_z=0;}

    /// Default destructor
    virtual ~SmearMCVertex(){}

    /** IMPLEMENT in SmearMCVertex.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SmearMCVertex.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SmearMCVertex.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetXOffset(double d) { _offset = d; }
    void FilterEvents(bool on) { _filter = on; }

  protected:

    bool   _filter;
    double _offset;
    double _time2cm;

    larutil::SpaceChargeMicroBooNE *_SCE;

    TH1D * _hist_x ;
    TH1D * _hist_y ;
    TH1D * _hist_z ;

    
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
