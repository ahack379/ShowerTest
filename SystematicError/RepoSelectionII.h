/**
 * \file RepoSelectionII.h
 *
 * \ingroup SystematicError
 * 
 * \brief Class def header for a class RepoSelectionII
 *
 * @author Code translated to larlite from SelectionII code in larsoft repository
 */

/** \addtogroup SystematicError

    @{*/

#ifndef LARLITE_REPOSELECTIONII_H
#define LARLITE_REPOSELECTIONII_H

#include "Analysis/ana_base.h"
#include "LArUtil/Geometry.h"

namespace larlite {
  /**
     \class RepoSelectionII
     User custom analysis class made by SHELL_USER_NAME
   */
  class RepoSelectionII : public ana_base{
  
  public:

    /// Default constructor
    RepoSelectionII(); //{ _name="RepoSelectionII"; _fout=0;}

    /// Default destructor
    virtual ~RepoSelectionII(){}

    /** IMPLEMENT in RepoSelectionII.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in RepoSelectionII.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in RepoSelectionII.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    bool inFV( double x, double y, double z) const; 

    double scaledEdx(double x, int plane, bool isdata) const;
     
    double random(double low, double high) ;

    void genNewValues() ;
    void SetDistToEdgeXYZ(double x, double y, double z){ 
        fDistToEdgeX = x;             
	fDistToEdgeY = y;
	fDistToEdgeZ = z;
	} 
    void SetPEThresh (double PE ) { fPEThresh = PE ; }
    void SetTrk2FlashDist (double t){ fTrk2FlashDist = t; }
    void SetMinTrk2VtxDist(double t){ fMinTrk2VtxDist = t;}
    void SetMinTrackLen        (double t){fMinTrackLen = t;} 
    void SetMaxCosineAngle     (double t){fMaxCosineAngle = t;} 
    void SetMaxCosy1stTrk      (double t){fMaxCosy1stTrk = t;} 
    void SetMinTrackLen2ndTrk  (double t){fMinTrackLen2ndTrk = t;} 
    void SetMaxCosySingle      (double t){fMaxCosySingle = t;} 
    void SetMinTrackLenSingle  (double t){fMinTrackLenSingle = t;} 
    void SetMindEdxRatioSingle (double t){fMindEdxRatioSingle = t;} 
    void SetMaxTrkLengthySingle(double t){fMaxTrkLengthySingle= t;} 
    void SetMinStartdEdx1stTrk (double t){fMinStartdEdx1stTrk = t;} 
    void SetMaxEnddEdx1stTrk   (double t){fMaxEnddEdx1stTrk = t;} 

  protected:

    std::string fTrackModuleLabel       ; 
    std::string fVertexModuleLabel      ; 
    std::string fOpFlashModuleLabel     ; 
    std::string fCalorimetryModuleLabel ; 
      
    double fDistToEdgeX             ;
    double fDistToEdgeY             ;
    double fDistToEdgeZ             ;
    double fBeamMin                 ; 
    double fBeamMax                 ;
    double fPEThresh                ;
    double fTrk2FlashDist           ;
    double fMinTrk2VtxDist          ;
    double fMinTrackLen             ;
    double fMaxCosineAngle          ;
    double fMaxCosy1stTrk           ;
    double fMinTrackLen2ndTrk       ;
    double fMaxCosySingle           ;
    double fMinTrackLenSingle       ;
    double fMindEdxRatioSingle      ;
    double fMaxTrkLengthySingle     ;
    double fMinStartdEdx1stTrk      ;
    double fMaxEnddEdx1stTrk        ;
    int    fDebug                   ;
    
    // Declare min and max parameters for random variable distribution
    double fMin_DistToEdgeX             ;
    double fMin_DistToEdgeY             ;
    double fMin_DistToEdgeZ             ;
    double fMin_BeamMin                 ; 
    double fMin_BeamMax                 ;
    double fMin_PEThresh                ;
    double fMin_Trk2FlashDist           ;
    double fMin_MinTrk2VtxDist          ;
    double fMin_MinTrackLen             ;
    double fMin_MaxCosineAngle          ;
    double fMin_MaxCosy1stTrk           ;
    double fMin_MinTrackLen2ndTrk       ;
    double fMin_MaxCosySingle           ;
    double fMin_MinTrackLenSingle       ;
    double fMin_MindEdxRatioSingle      ;
    double fMin_MaxTrkLengthySingle     ;
    double fMin_MinStartdEdx1stTrk      ;
    double fMin_MaxEnddEdx1stTrk        ;

    // Declare min and max parameters for random variable distribution
    double fMax_DistToEdgeX             ;
    double fMax_DistToEdgeY             ;
    double fMax_DistToEdgeZ             ;
    double fMax_BeamMin                 ; 
    double fMax_BeamMax                 ;
    double fMax_PEThresh                ;
    double fMax_Trk2FlashDist           ;
    double fMax_MinTrk2VtxDist          ;
    double fMax_MinTrackLen             ;
    double fMax_MaxCosineAngle          ;
    double fMax_MaxCosy1stTrk           ;
    double fMax_MinTrackLen2ndTrk       ;
    double fMax_MaxCosySingle           ;
    double fMax_MinTrackLenSingle       ;
    double fMax_MindEdxRatioSingle      ;
    double fMax_MaxTrkLengthySingle     ;
    double fMax_MinStartdEdx1stTrk      ;
    double fMax_MaxEnddEdx1stTrk        ;

    const larutil::Geometry * fGeometry ;

    /////////// My extra variables
    int _n_it_per_event ;
    
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
