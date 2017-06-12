/**
 * \file AcceptanceStudy.h
 *
 * \ingroup CalcEfficiency 
 * 
 * \brief Class def header for a class AcceptanceStudy
 *
 * @author ah673
 */

/** \addtogroup CalcEfficiency 

    @{*/

#ifndef LARLITE_ACCEPTANCESTUDY_H
#define LARLITE_ACCEPTANCESTUDY_H

#include "Analysis/ana_base.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/SpaceChargeMicroBooNE.h"
#include "GeoAlgo/GeoAlgo.h"

namespace larlite {
  /**
     \class AcceptanceStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class AcceptanceStudy : public ana_base{
  
  public:

    /// Default constructor
    AcceptanceStudy(){ _name="AcceptanceStudy"; _fout=0; }

    /// Default destructor
    virtual ~AcceptanceStudy(){}

    /** IMPLEMENT in AcceptanceStudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in AcceptanceStudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in AcceptanceStudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  std::multimap<float,float> _map;
  int _event ; 

  int _infv_ccpi0 ; 
  int _thresh ; // Study how many ccpi0 below threshold
  int _dalitz ;
  int _out_of_vol;
  
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
