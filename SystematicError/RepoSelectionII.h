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
#include "GeoAlgo/GeoAlgo.h"
#include "LArUtil/GeometryHelper.h"

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

    void genNewSelIIValues() ;

    void genNewHitRatioValues() ;

    void genNewShowerPairValues() ;

    /// Set SelectionII Variables
    void SetBounds_DistToEdgeXYZ(double minx, double miny, double minz,
                                 double maxx, double maxy, double maxz){ 
        fMin_DistToEdgeX = minx;             
	fMin_DistToEdgeY = miny;
	fMin_DistToEdgeZ = minz;
        fMax_DistToEdgeX = maxx;             
	fMax_DistToEdgeY = maxy;
	fMax_DistToEdgeZ = maxz;
	} 
    void SetBounds_PEThresh           (double min, double max){fMin_PEThresh           = min; fMax_PEThresh           = max;}  
    void SetBounds_Trk2FlashDist      (double min, double max){fMin_Trk2FlashDist      = min; fMax_Trk2FlashDist      = max;}
    void SetBounds_MinTrk2VtxDist     (double min, double max){fMin_MinTrk2VtxDist     = min; fMax_MinTrk2VtxDist     = max;}
    void SetBounds_MinTrackLen        (double min, double max){fMin_MinTrackLen        = min; fMax_MinTrackLen        = max;}
    void SetBounds_MaxCosineAngle     (double min, double max){fMin_MaxCosineAngle     = min; fMax_MaxCosineAngle     = max;}
    void SetBounds_MaxCosy1stTrk      (double min, double max){fMin_MaxCosy1stTrk      = min; fMax_MaxCosy1stTrk      = max;}
    void SetBounds_MinTrackLen2ndTrk  (double min, double max){fMin_MinTrackLen2ndTrk  = min; fMax_MinTrackLen2ndTrk  = max;}
    void SetBounds_MaxCosySingle      (double min, double max){fMin_MaxCosySingle      = min; fMax_MaxCosySingle      = max;}
    void SetBounds_MinTrackLenSingle  (double min, double max){fMin_MinTrackLenSingle  = min; fMax_MinTrackLenSingle  = max;}
    void SetBounds_MindEdxRatioSingle (double min, double max){fMin_MindEdxRatioSingle = min; fMax_MindEdxRatioSingle = max;}
    void SetBounds_MaxTrkLengthySingle(double min, double max){fMin_MaxTrkLengthySingle= min; fMax_MaxTrkLengthySingle= max;}
    void SetBounds_MinStartdEdx1stTrk (double min, double max){fMin_MinStartdEdx1stTrk = min; fMax_MinStartdEdx1stTrk = max;}
    void SetBounds_MaxEnddEdx1stTrk   (double min, double max){fMin_MaxEnddEdx1stTrk   = min; fMax_MaxEnddEdx1stTrk   = max;}

    /// Set HitRatio Variables
    void SetBounds_Radius        (float min, float max){fMin_Radius          = min; fMax_Radius          = max;} 
    void SetBounds_HitRatio      (float min, float max){fMin_HitRatio        = min; fMax_HitRatio        = max;}
    void SetBounds_MinHitsRequired (int min  , int max  ){fMin_MinHitsRequired = min; fMax_MinHitsRequired = max;}

    /// Set Shower Variables
    void SetBounds_OpeningAngle  (float min, float max) { fMin_MinOpeningAngle = min; fMax_MinOpeningAngle = max;} 
    void SetBounds_MaxIP         (float min, float max) { fMin_MaxIP           = min; fMax_MaxIP           = max;}
    void SetBounds_MaxRadLength  (float min, float max) { fMin_MaxRadLength    = min; fMax_MaxRadLength    = max;}

    /////////////////////////////////////
    // Other Variables
    ////////////////////////////////////
    void SetIters ( int it ) { _n_it_per_event = it; }
    
  protected:
    
    /////////////////////////////
    // Selection II variables  //
    /////////////////////////////
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

    /////////////////////////////
    // HitRatio variables      //
    /////////////////////////////
    float fRadius ;
    float fHitRatio ;
    int   fMinHitsRequired ;

    // Min and max parameters for random generation
    float fMin_Radius ;
    float fMin_HitRatio ;
    float fMin_MinHitsRequired ;
    float fMax_Radius ;
    float fMax_HitRatio ;
    float fMax_MinHitsRequired ;

    /////////////////////////////
    // Showerreco variables    //
    /////////////////////////////
    float fMinOpeningAngle ;
    float fMaxIP ;
    float fMaxRadLength ;

    // Min and max parameters for random generation
    float fMin_MinOpeningAngle ;
    float fMin_MaxIP ;
    float fMin_MaxRadLength ;

    float fMax_MinOpeningAngle ;
    float fMax_MaxIP ;
    float fMax_MaxRadLength ;

    ::geoalgo::GeoAlgo _geoAlgo ;
    const larutil::GeometryHelper * fGeomH ;


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
