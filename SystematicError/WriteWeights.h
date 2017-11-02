/**
 * \file WriteWeights.h
 *
 * \ingroup Playground
 * 
 * \brief Class def header for a class WriteWeights
 *
 * @author davidc1
 */

/** \addtogroup Playground

    @{*/

#ifndef LARLITE_WRITEWEIGHTS_H
#define LARLITE_WRITEWEIGHTS_H

#include "Analysis/ana_base.h"
#include <fstream>

namespace larlite {
  /**
     \class WriteWeights
     User custom analysis class made by SHELL_USER_NAME
   */
  class WriteWeights : public ana_base{
  
  public:

    /// Default constructor
    WriteWeights(){ _name="WriteWeights"; _fout=0; _event_producer = "genieeventweight"; }

    /// Default destructor
    virtual ~WriteWeights(){}

    /** IMPLEMENT in WriteWeights.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in WriteWeights.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in WriteWeights.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetWeightProducer( std::string prod ) { _event_producer = prod ; }

  protected:

    std::ofstream event_file;
    int _event ;

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
