/**
 * \file CalcCharge.h
 *
 * \ingroup SumCharge
 * 
 * \brief Class def header for a class CalcCharge
 *
 * @author ah673
 */

/** \addtogroup SumCharge

    @{*/

#ifndef LARLITE_CALCCHARGE_H
#define LARLITE_CALCCHARGE_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class CalcCharge
     User custom analysis class made by SHELL_USER_NAME
   */
  class CalcCharge : public ana_base{
  
  public:

    /// Default constructor
    CalcCharge(){ _name="CalcCharge"; _fout=0; _tree=nullptr ;}

    /// Default destructor
    virtual ~CalcCharge(){}

    /** IMPLEMENT in CalcCharge.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CalcCharge.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CalcCharge.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  TTree * _tree ;
  float _sum_charge00 ;
  float _sum_charge02 ;
    
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
