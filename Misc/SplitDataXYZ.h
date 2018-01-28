/**
 * \file SplitDataXYZ.h
 *
 * \ingroup Misc
 * 
 * \brief Class def header for a class SplitDataXYZ
 *
 * @author ah673
 */

/** \addtogroup Misc

    @{*/

#ifndef LARLITE_SPLITDATAXYZ_H
#define LARLITE_SPLITDATAXYZ_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class SplitDataXYZ
     User custom analysis class made by SHELL_USER_NAME
   */
  class SplitDataXYZ : public ana_base{
  
  public:

    /// Default constructor
    SplitDataXYZ(){ _name="SplitDataXYZ"; _fout=0; _tree=0;}

    /// Default destructor
    virtual ~SplitDataXYZ(){}

    /** IMPLEMENT in SplitDataXYZ.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SplitDataXYZ.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SplitDataXYZ.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:

    TTree * _tree ;
  
    int _x_pos;
    int _y_pos;
    int _z_pos;
    int _x_neg;
    int _y_neg;
    int _z_neg;

    float _x_pos_POT;
    float _y_pos_POT;
    float _z_pos_POT;
    float _x_neg_POT;
    float _y_neg_POT;
    float _z_neg_POT;


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
