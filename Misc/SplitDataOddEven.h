/**
 * \file SplitDataOddEven.h
 *
 * \ingroup Misc
 * 
 * \brief Class def header for a class SplitDataOddEven
 *
 * @author ah673
 */

/** \addtogroup Misc

    @{*/

#ifndef LARLITE_SPLITDATAODDEVEN_H
#define LARLITE_SPLITDATAODDEVEN_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class SplitDataOddEven
     User custom analysis class made by SHELL_USER_NAME
   */
  class SplitDataOddEven : public ana_base{
  
  public:

    /// Default constructor
    SplitDataOddEven(){ _name="SplitDataOddEven"; _fout=0; _tree=0;}

    /// Default destructor
    virtual ~SplitDataOddEven(){}

    /** IMPLEMENT in SplitDataOddEven.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SplitDataOddEven.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SplitDataOddEven.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TTree * _tree ;
  
    int _x_odd;
    int _x_even;

    float _x_odd_POT;
    float _x_even_POT;


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
