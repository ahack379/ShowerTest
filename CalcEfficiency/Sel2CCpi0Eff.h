/**
 * \file Sel2CCpi0Eff.h
 *
 * \ingroup CalcEfficiency
 * 
 * \brief Class def header for a class Sel2CCpi0Eff
 *
 * @author ah673
 */

/** \addtogroup CalcEfficiency

    @{*/

#ifndef LARLITE_SEL2CCPI0EFF_H
#define LARLITE_SEL2CCPI0EFF_H

#include "Analysis/ana_base.h"
#include <TTree.h>

namespace larlite {
  /**
     \class Sel2CCpi0Eff
     User custom analysis class made by SHELL_USER_NAME
   */
  class Sel2CCpi0Eff : public ana_base{
  
  public:

    /// Default constructor
    Sel2CCpi0Eff(){ _name="Sel2CCpi0Eff"; _fout=0; _CCNC=false; _cut_tree=nullptr;}

    /// Default destructor
    virtual ~Sel2CCpi0Eff(){}

    /** IMPLEMENT in Sel2CCpi0Eff.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in Sel2CCpi0Eff.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void GetBothCCNC( bool getem){ _CCNC = getem ; }

  protected:

  int _events ;
  int _signal ;

  bool _CCNC ;

  std::vector<int> _event_list ;

  TTree * _cut_tree ;
  float _energy ;
    
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
