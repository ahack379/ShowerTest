/**
 * \file ShowerResolution.h
 *
 * \ingroup EnergyStudy
 * 
 * \brief Class def header for a class ShowerResolution
 *
 * @author ah673
 */

/** \addtogroup EnergyStudy

    @{*/

#ifndef LARLITE_SHOWERRESOLUTION_H
#define LARLITE_SHOWERRESOLUTION_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class ShowerResolution
     User custom analysis class made by SHELL_USER_NAME
   */
  class ShowerResolution : public ana_base{
  
  public:

    /// Default constructor
    ShowerResolution(){ _name="ShowerResolution"; _fout=0; _tree = 0 ;}

    /// Default destructor
    virtual ~ShowerResolution(){}

    /** IMPLEMENT in ShowerResolution.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ShowerResolution.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ShowerResolution.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

  TTree * _tree; 
  float _e_mcc;
  float _e_true;

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
