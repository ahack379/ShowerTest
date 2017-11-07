/**
 * \file SaveWeights.h
 *
 * \ingroup Playground
 * 
 * \brief Class def header for a class SaveWeights
 *
 * @author davidc1
 */

/** \addtogroup Playground

    @{*/

#ifndef LARLITE_SAVEWEIGHTS_H
#define LARLITE_SAVEWEIGHTS_H

#include "Analysis/ana_base.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>

namespace larlite {
  /**
     \class SaveWeights
     User custom analysis class made by SHELL_USER_NAME
   */
  class SaveWeights : public ana_base{
  
  public:

    /// Default constructor
    SaveWeights(){ _name="SaveWeights"; _fout=0; _wgtmap.clear(); _event_producer = "genieeventweight"; }

    /// Default destructor
    virtual ~SaveWeights(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void SetWeightProducer( std::string prod ) { _event_producer = prod ; }

    //void addWeight(const int& run, const int & subrun, const int& event, double *) ; //std::vector<double> ); 

  protected:
    
    std::ifstream _file;

    std::map< std::pair<int,int> , std::vector<double> > _wgtmap;

    std::vector<std::string> _func_v ;

    std::string _event_producer ;
    
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
