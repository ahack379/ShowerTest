/**
 * \file XSecError.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class XSecError
 *
 * @author ah673
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_XSECERROR_H
#define LARLITE_XSECERROR_H

#include "Analysis/ana_base.h"
#include "TH1F.h"

namespace larlite {
  /**
     \class XSecError
     User custom analysis class made by SHELL_USER_NAME
   */
  class XSecError : public ana_base{
  
  public:

    /// Default constructor
    XSecError(){ _name="XSecError"; _fout=0; _tree=nullptr; _hist=0; _hist2 = 0; _hist3=0; _histw = 0;}

    /// Default destructor
    virtual ~XSecError(){}

    /** IMPLEMENT in XSecError.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in XSecError.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in XSecError.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  TTree * _tree;
  TH1F * _hist;
  TH1F * _hist2;
  TH1F * _hist3;
  TH1F * _histw;
  float _test; 
  float _test_weight ;
    
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
