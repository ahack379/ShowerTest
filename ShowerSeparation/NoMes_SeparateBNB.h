/**
 * \file NoMes_SeparateBNB.h
 *
 * \ingroup HitDensity
 * 
 * \brief Class def header for a class NoMes_SeparateBNB
 *
 * @author ah673
 */

/** \addtogroup HitDensity

    @{*/

#ifndef LARLITE_NOMES_SEPARATEBNB_H
#define LARLITE_NOMES_SEPARATEBNB_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class NoMes_SeparateBNB
     User custom analysis class made by SHELL_USER_NAME
   */
  class NoMes_SeparateBNB : public ana_base{
  
  public:

    /// Default constructor
    NoMes_SeparateBNB(){ _name="NoMes_SeparateBNB"; _fout=0; _get_pi0s = false ;}

    /// Default destructor
    virtual ~NoMes_SeparateBNB(){}

    /** IMPLEMENT in NoMes_SeparateBNB.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in NoMes_SeparateBNB.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in NoMes_SeparateBNB.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void GetPi0s(bool get){ _get_pi0s = get ; }

  protected:

  bool _get_pi0s;
  int _event ;
  int _signal;
    
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
