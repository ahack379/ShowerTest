/**
 * \file MakeVtx.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class MakeVtx
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_MAKEVTX_H
#define LARLITE_MAKEVTX_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class MakeVtx
     User custom analysis class made by SHELL_USER_NAME
   */
  class MakeVtx : public ana_base{
  
  public:

    /// Default constructor
    MakeVtx(){ _name="MakeVtx"; _fout=0; _pi0=1; _pi0_offset=0;_mu_offset = 0;}

    /// Default destructor
    virtual ~MakeVtx(){}

    /** IMPLEMENT in MakeVtx.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MakeVtx.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MakeVtx.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetPi0( bool usingPi0 ) { _pi0 = usingPi0 ; }

    void SetPi0Offset(float offset) { _pi0_offset = offset ; } 
    void SetMuOffset(float offset) { _mu_offset = offset ; } 


  protected:

  int _id; 

  bool _pi0 ;
  float _pi0_offset ;
  float _mu_offset ;
    
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
