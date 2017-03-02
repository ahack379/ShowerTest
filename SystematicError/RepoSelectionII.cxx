#ifndef LARLITE_REPOSELECTIONII_CXX
#define LARLITE_REPOSELECTIONII_CXX

#include "RepoSelectionII.h"
#include "DataFormat/track.h"
#include "DataFormat/opflash.h"
#include "DataFormat/vertex.h"
#include "DataFormat/calorimetry.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/hit.h"
#include "DataFormat/shower.h"
#include "DataFormat/potsummary.h"
#include "DataFormat/mctruth.h"

// HitRatio
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

namespace larlite {

  RepoSelectionII::RepoSelectionII() {
   
    fTrackModuleLabel        = "pandoraNu"; 
    fVertexModuleLabel       = "pandoraNu";
    fOpFlashModuleLabel      = "simpleFlashBeam"; 
    fCalorimetryModuleLabel  = "pandoraNucalo";
   
    ///////////////////////////////////////////////////
    // Default values in larsoft updated SelectionII //
    ///////////////////////////////////////////////////
    fDistToEdgeX             =   20.;
    fDistToEdgeY             =   20.;
    fDistToEdgeZ             =   10.;
    fBeamMin                 =  3.55 ; //these values are for mcc7 beam window; may need to adjust 3.3;
    fBeamMax                 =  5.15 ; //these values are for mcc7 beam window; may need to adjust 4.9;
    fPEThresh                =   50.;
    fTrk2FlashDist           =   70.;
    fMinTrk2VtxDist          =    3.;
    fMinTrackLen             =   15.;
    fMaxCosineAngle          =   0.9;
    fMaxCosy1stTrk           =   0.6;
    fMinTrackLen2ndTrk       =   30.;
    fMaxCosySingle           =   0.7;
    fMinTrackLenSingle       =   40.;
    fMindEdxRatioSingle      =   1.5;
    fMaxTrkLengthySingle     =   25.;
    fMinStartdEdx1stTrk      =   2.5;
    fMaxEnddEdx1stTrk        =   4.0;
    fDebug                   =     0;
    
    // Min on each parameter
    fMin_DistToEdgeX             =   5.;
    fMin_DistToEdgeY             =   5.;
    fMin_DistToEdgeZ             =   0.;
    fMin_PEThresh                =   40.;
    fMin_Trk2FlashDist           =   55.;
    fMin_MinTrk2VtxDist          =    0.;
    fMin_MinTrackLen             =   5.;
    fMin_MaxCosineAngle          =   0.85;
    fMin_MaxCosy1stTrk           =   0.5;
    fMin_MinTrackLen2ndTrk       =   20.;
    fMin_MaxCosySingle           =   0.5;
    fMin_MinTrackLenSingle       =   30.;
    fMin_MindEdxRatioSingle      =   1.;
    fMin_MaxTrkLengthySingle     =   15.;
    fMin_MinStartdEdx1stTrk      =   2.1;
    fMin_MaxEnddEdx1stTrk        =   1.0;

    // Max on each parameter 
    fMax_DistToEdgeX             =   38.;
    fMax_DistToEdgeY             =   38.;
    fMax_DistToEdgeZ             =   20.;
    fMax_PEThresh                =   60.;
    fMax_Trk2FlashDist           =   85.;
    fMax_MinTrk2VtxDist          =    7.;
    fMax_MinTrackLen             =   50.;
    fMax_MaxCosineAngle          =   1.;
    fMax_MaxCosy1stTrk           =   0.85;
    fMax_MinTrackLen2ndTrk       =   50.;
    fMax_MaxCosySingle           =   0.7;
    fMax_MinTrackLenSingle       =   40.;
    fMax_MindEdxRatioSingle      =   2. ;
    fMax_MaxTrkLengthySingle     =   45.;
    fMax_MinStartdEdx1stTrk      =   3.5;
    fMax_MaxEnddEdx1stTrk        =   10.0;

    ////////////////////////////////////////////
    // Default values for HitRatio            //
    ////////////////////////////////////////////
    fRadius                      = 60. ;
    fHitRatio                    = 0.24 ;
    fMinHitsRequired             = 15 ;

    // Min and Max on each parameter 
    fMin_Radius                  = 40. ;
    fMin_HitRatio                = 0.15 ;
    fMin_MinHitsRequired         = 10 ;

    fMax_Radius                  = 75. ;
    fMax_HitRatio                = 0.38 ;
    fMax_MinHitsRequired         = 50 ;

    ////////////////////////////////////////////
    // Default values for Showerreco          //
    ////////////////////////////////////////////
    fMinOpeningAngle             = 0.35 ;
    fMaxIP                       = 4. ;
    fMaxRadLength                = 62 ;

    // Min and Max on each parameter 
    fMin_MinOpeningAngle         = 0.3 ;
    fMin_MaxIP                   = 2. ;
    fMin_MaxRadLength            = 40 ;

    fMax_MinOpeningAngle         = 0.65 ;
    fMax_MaxIP                   = 8. ;
    fMax_MaxRadLength            = 80 ;
    
    _name                    = "RepoSelectionII";
    _fout                    = 0;
    fGeometry = nullptr;
    fGeomH = nullptr;
    _n_it_per_event = 100;
    _it_tree = nullptr;

  }
  
  bool RepoSelectionII::initialize(){

    fGeometry = larutil::Geometry::GetME(); 
    fGeomH = larutil::GeometryHelper::GetME(); 

    _n_true_signal = 0;
    _tot_pot = 0. ;

    if( !_it_tree ){
       _it_tree = new TTree("it_tree","it_tree");
       _it_tree->Branch("tagged_v",&_tagged_v,"std::vector<int>");
       _it_tree->Branch("bkgd_v",&_tagged_v,"std::vector<int>");
       _it_tree->Branch("eff_v",&_tagged_v,"std::vector<int>");
       _it_tree->Branch("xsec_v",&_tagged_v,"std::vector<int>");
       }

    _tagged_v.resize(_n_it_per_event,0) ;
    _bkgd_v.resize(_n_it_per_event,0) ;
    _eff_v.resize(_n_it_per_event,0) ;
    _xsec_v.resize(_n_it_per_event,0) ;
 
    return true;
   }

  void RepoSelectionII::genNewSelIIValues(){

    //fDistToEdgeX = random(fMin_DistToEdgeX,fMax_DistToEdgeX);
    //fDistToEdgeY = random(fMin_DistToEdgeY,fMax_DistToEdgeY);
    //fDistToEdgeZ = random(fMin_DistToEdgeZ,fMax_DistToEdgeZ);
    fPEThresh    = random(fMin_PEThresh,fMax_PEThresh);
    fTrk2FlashDist  = random(fMin_Trk2FlashDist,fMax_Trk2FlashDist);
    fMinTrk2VtxDist = random(fMin_MinTrk2VtxDist,fMax_MinTrk2VtxDist);
    fMinTrackLen    = random(fMin_MinTrackLen,fMax_MinTrackLen);
    fMaxCosineAngle = random(fMin_MaxCosineAngle,fMax_MaxCosineAngle); 
    fMaxCosy1stTrk  = random(fMin_MaxCosy1stTrk,fMax_MaxCosy1stTrk);
    fMinTrackLen2ndTrk  = random(fMin_MinTrackLen2ndTrk,fMax_MinTrackLen2ndTrk);
    fMaxCosySingle      = random(fMin_MaxCosySingle,fMax_MaxCosySingle) ;
    fMinTrackLenSingle  = random(fMin_MinTrackLenSingle,fMax_MinTrackLenSingle);
    fMindEdxRatioSingle = random(fMin_MindEdxRatioSingle,fMax_MindEdxRatioSingle);
    fMaxTrkLengthySingle = random(fMin_MaxTrkLengthySingle,fMax_MaxTrkLengthySingle);
    fMinStartdEdx1stTrk  = random(fMin_MinStartdEdx1stTrk,fMax_MinStartdEdx1stTrk);
    fMaxEnddEdx1stTrk    = random(fMin_MaxEnddEdx1stTrk,fMax_MaxEnddEdx1stTrk);

   }

  void RepoSelectionII::genNewHitRatioValues(){

    fRadius = random(fMin_Radius,fMax_Radius) ;
    fHitRatio = random(fMin_HitRatio,fMax_HitRatio) ;
    fMinHitsRequired = random(fMin_MinHitsRequired,fMax_MinHitsRequired) ;

     }

  void RepoSelectionII::genNewShowerPairValues(){

    fMinOpeningAngle = random(fMin_MinOpeningAngle,fMin_MinOpeningAngle) ;
    fMaxIP = random(fMin_MaxIP,fMax_MaxIP) ;
    fMaxRadLength = random(fMin_MaxRadLength,fMax_MaxRadLength) ;

     }
 
  bool RepoSelectionII::analyze(storage_manager* storage) {

    auto flashlist = storage->get_data<event_opflash>(fOpFlashModuleLabel);
    auto tracklist = storage->get_data<event_track>(fTrackModuleLabel);
    auto vtxlist = storage->get_data<event_vertex>(fVertexModuleLabel);
    auto calolist = storage->get_data<event_calorimetry>(fCalorimetryModuleLabel);
    auto ev_hit_g = storage->get_data<event_hit>("gaushit");
    auto ev_s = storage->get_data<event_shower>("showerreco");

    if( !flashlist || flashlist->size() == 0 ){ 
       std::cout<<"No flashes..."<<std::endl ;
       return false;
       }

    if( !tracklist || tracklist->size() == 0 ){ 
       std::cout<<"No tracks..."<<std::endl ;
       return false;
       }

    if( !vtxlist || vtxlist->size() == 0 ){ 
       std::cout<<"No vertices..."<<std::endl ;
       return false;
       }
  
    if( !calolist || calolist->size() == 0 ){ 
       std::cout<<"No calo..."<<std::endl ;
       return false;
       }

    if ( !ev_hit_g || !ev_hit_g->size() ) {
       std::cout<<"No hits, returning..."<<std::endl ; 
       return false; 
    }

    if( !ev_s || !ev_s->size() || ev_s->size() < 2 ){
      std::cout<<"Not enough reco'd showers..." <<std::endl;
      return false;
     }
  
    //art::FindMany<anab::Calorimetry>  fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    auto ev_ass = storage->get_data<larlite::event_ass>(fCalorimetryModuleLabel);
    auto const& ass_calo_v = ev_ass->association(tracklist->id(),calolist->id());

    if( !ev_ass || ev_ass->size() == 0 ){ 
       std::cout<<"No calo..."<<std::endl ;
       return false;
       }
  
    if( ass_calo_v.size() == 0 ){ 
       std::cout<<"No calo ass..."<<std::endl ;
       return false;
       }

    auto ev_pot = storage->get_subrundata<potsummary>("generator"); 
    auto ev_mctruth = storage->get_data<event_mctruth>("generator"); 

    if( !ev_pot ){ 
       std::cout<<"No POT..."<<std::endl ;
       return false;
       }

    if(!ev_mctruth || !ev_mctruth->size() ){ 
       std::cout<<"No Truth..."<<std::endl ;
       return false;
    }

    if( storage->subrun_id() != storage->last_subrun_id() )
      _tot_pot += ev_pot->totgoodpot ;

    auto nu  = ev_mctruth->at(0).GetNeutrino();
    double FV_test_xyz[3] = {0.};
    auto traj = nu.Nu().Trajectory();
    FV_test_xyz[0] = traj.at(traj.size() - 1).X();
    FV_test_xyz[1] = traj.at(traj.size() - 1).Y();
    FV_test_xyz[2] = traj.at(traj.size() - 1).Z();

    auto parts = ev_mctruth->at(0).GetParticles();

    int n_pi0 = 0;
    int n_mu = 0;
    bool isSignal = false ;
    
    for ( auto const & p : parts ){
    
      if( p.StatusCode() == 1 && p.PdgCode() == 111 )
        n_pi0 += 1;

      if( p.StatusCode() == 1 && p.PdgCode() == 13 )
        n_mu += 1;
        }

    if( (FV_test_xyz[0] > 0 + fDistToEdgeX || FV_test_xyz[0] < 256.35 - fDistToEdgeX || 
          FV_test_xyz[1] > -116.5 + fDistToEdgeY ||  FV_test_xyz[1] < 116.5 - fDistToEdgeY || 
          FV_test_xyz[2] > 0 + fDistToEdgeZ || FV_test_xyz[2] < 1036.8 - fDistToEdgeZ) && 
          n_mu == 1 && n_pi0 == 1 ){
        isSignal = true; 
        _n_true_signal ++ ;
        }


    for( int evt_i = 0; evt_i < _n_it_per_event; evt_i++){
       
        if ( evt_i != 0 ) genNewSelIIValues();

        //check the flash info
        double FlashPEmax=0;
        int NuFlashID=-1;

        for (size_t i = 0; i<flashlist->size(); ++i){
          auto f_i = flashlist->at(i) ;
          if (f_i.TotalPE()>fPEThresh && f_i.Time()>fBeamMin && f_i.Time()<fBeamMax){
            if (f_i.TotalPE()>FlashPEmax){
              FlashPEmax = f_i.TotalPE();
              NuFlashID = i;
            }
          }
        }

        //Did not find the desired flash, return
        if (NuFlashID == -1)
          return false;

        //Save basic track information in vectors
        std::vector<double> trkstartx(tracklist->size());
        std::vector<double> trkstarty(tracklist->size());
        std::vector<double> trkstartz(tracklist->size());
        std::vector<double> trkendx(tracklist->size());
        std::vector<double> trkendy(tracklist->size());
        std::vector<double> trkendz(tracklist->size());
        std::vector<double> trkstartdcosx(tracklist->size());
        std::vector<double> trkstartdcosy(tracklist->size());
        std::vector<double> trkstartdcosz(tracklist->size());
        std::vector<double> trkenddcosx(tracklist->size());
        std::vector<double> trkenddcosy(tracklist->size());
        std::vector<double> trkenddcosz(tracklist->size());
        std::vector<double> trklen(tracklist->size());
        double larStart[3];
        double larEnd[3];
        std::vector<double> trackStart;
        std::vector<double> trackEnd;
        for (size_t i = 0; i<tracklist->size(); ++i){
          auto t_i = tracklist->at(i) ;

          trackStart.clear();
          trackEnd.clear();
          memset(larStart, 0, 3);
          memset(larEnd, 0, 3);
          t_i.Extent(trackStart,trackEnd);
          t_i.Direction(larStart,larEnd);
          trkstartx[i]      = trackStart[0];
          trkstarty[i]      = trackStart[1];
          trkstartz[i]      = trackStart[2];
          trkendx[i]        = trackEnd[0];
          trkendy[i]        = trackEnd[1];
          trkendz[i]        = trackEnd[2];
          trkstartdcosx[i]  = larStart[0];
          trkstartdcosy[i]  = larStart[1];
          trkstartdcosz[i]  = larStart[2];
          trkenddcosx[i]    = larEnd[0];
          trkenddcosy[i]    = larEnd[1];
          trkenddcosz[i]    = larEnd[2];
          trklen[i]         = t_i.Length();
        }

        //Match each track with the selected flash
        std::vector<bool> trackflashmatch(tracklist->size());
        bool foundtrackflashmatch = false;

        //double TaggedFlashYCenter = flashlist[NuFlashID]->YCenter();
        double TaggedFlashZCenter = flashlist->at(NuFlashID).ZCenter();
        for (size_t i = 0; i<tracklist->size(); ++i){
          double FlashTrackDis = 1e10;
          if ((trkstartz[i]<TaggedFlashZCenter && trkendz[i]>TaggedFlashZCenter)||
              (trkstartz[i]>TaggedFlashZCenter && trkendz[i]<TaggedFlashZCenter)){
            FlashTrackDis = 0;
          }
          else{
            FlashTrackDis = std::min(std::abs(trkstartz[i] - TaggedFlashZCenter),
                                     std::abs(trkendz[i] - TaggedFlashZCenter));
          }
          if (FlashTrackDis<fTrk2FlashDist){
            trackflashmatch[i] = true;
            foundtrackflashmatch = true;
          }
          else{
            trackflashmatch[i] = false;
          }
        }
        if (!foundtrackflashmatch) {
          if (fDebug) std::cout<<"Did not find any tracks matching flash."<<std::endl;
          return false;
        }

        //Match tracks with vertices
        if (!vtxlist->size()) {
          if (fDebug) std::cout<<"No vertex found"<<std::endl;
          return false;
        }

        std::vector<std::vector<int>> trkindex(vtxlist->size());
        std::vector<std::vector<bool>> fliptrack(vtxlist->size());
        for (size_t i = 0; i<vtxlist->size(); ++i){
          double xyz[3];
          vtxlist->at(i).XYZ(xyz);
          if (!inFV(xyz[0], xyz[1], xyz[2])) continue;
            for (size_t j = 0; j<tracklist->size(); ++j){
              double vtxtrkStartDis = sqrt(pow(trkstartx[j]-xyz[0],2)+
                                           pow(trkstarty[j]-xyz[1],2)+
                                           pow(trkstartz[j]-xyz[2],2));
              double vtxtrkEndDis = sqrt(pow(trkendx[j]-xyz[0],2)+
                                         pow(trkendy[j]-xyz[1],2)+
                                         pow(trkendz[j]-xyz[2],2));
              double vtxtrkDis = std::min(vtxtrkStartDis, vtxtrkEndDis);
              if (vtxtrkDis<fMinTrk2VtxDist){
                trkindex[i].push_back(j);
                if (vtxtrkEndDis<vtxtrkStartDis){
                  fliptrack[i].push_back(true);
                }
                else{
                  fliptrack[i].push_back(false);
                }
              }
            }//Loope over all tracks
        }//Loop over all vertices

        //calculate average dE/dx near the track start and track end
        std::vector<double> trkStartdEdx(tracklist->size());
        std::vector<double> trkEnddEdx(tracklist->size());
        for (size_t i = 0; i<tracklist->size(); ++i){
          //if (fmcal.isValid()){
            //std::vector<const anab::Calorimetry*> calos = fmcal.at(i);
            // calos here is Calo ID vector 
            auto calos = ass_calo_v.at(i);
            int icalo = -1;
            int totalnhits = 0;
            for (size_t j = 0; j<calos.size(); ++j){
              auto calo = calolist->at(j) ; 
              if (int(calo.dEdx().size())>totalnhits){
                icalo = j;
                totalnhits = calo.dEdx().size();
              }
            //}

            double sumdEdxStart=0;
            double sumdEdxEnd=0;

            int MaxHits=0;
            if(totalnhits>=20){
              MaxHits=10;
            }
            else if(totalnhits>0){
              MaxHits=totalnhits/2;
            }
            for(int ihit=0;ihit<MaxHits;ihit++){
              sumdEdxStart += calo.dEdx()[ihit]*scaledEdx(calo.XYZ()[ihit].X(), calo.PlaneID().Plane, 0) ; //evt.isRealData());
              sumdEdxEnd += calo.dEdx()[totalnhits-ihit-1]*scaledEdx(calo.XYZ()[totalnhits-ihit-1].X(), calo.PlaneID().Plane, 0); //evt.isRealData());
            }
           if (icalo!=-1&&
                sqrt(pow(calo.XYZ()[0].X()-trkstartx[i],2)+
                     pow(calo.XYZ()[0].Y()-trkstarty[i],2)+
                     pow(calo.XYZ()[0].Z()-trkstartz[i],2))>
                sqrt(pow(calo.XYZ()[0].X()-trkendx[i],2)+
                     pow(calo.XYZ()[0].Y()-trkendy[i],2)+
                     pow(calo.XYZ()[0].Z()-trkendz[i],2))){
              std::swap(sumdEdxEnd, sumdEdxStart);
            }
            if(MaxHits>0 && sumdEdxEnd>0 && sumdEdxStart>0){
              trkStartdEdx[i] = sumdEdxStart/MaxHits;
              trkEnddEdx[i] = sumdEdxEnd/MaxHits;
            }
            if (fDebug) std::cout<<"Trkid = "<<i<<" MaxHits = "<<MaxHits<<" sumdEdxStart "<<sumdEdxStart<<" sumdEdxEnd "<<sumdEdxEnd<<" icalo "<<icalo<<std::endl;
          }
        }//Loop over tracks

        //Examine tracks around each vertex and select neutrino candidates
        std::vector<bool> nuvtx(vtxlist->size());
        std::vector<double> cosangle(vtxlist->size()); //angle between two longest tracks
        std::vector<double> dcosylong(vtxlist->size()); //dcosy of the longest track
        std::vector<double> trklen2nd(vtxlist->size()); //track length of the second longest track
        for (size_t i = 0; i<vtxlist->size(); ++i){

          //Check if there are tracks associated with this vertex
          if (!trkindex[i].size()) {
              //no tracks associated with vertex
            nuvtx[i] = false;
            if (fDebug) std::cout<<"ivtx = "<<i<<" no tracks associated with this vertex."<<std::endl;
            continue;
          }

          //Check if at least one track matches the flash
          bool flashmatch = false;
          for (size_t j = 0; j<trkindex[i].size(); ++j){
            if (trackflashmatch[trkindex[i][j]]){
              flashmatch = true;
            }
          }
          if (!flashmatch){
            //no tracks matched to the flash around this vertex
              nuvtx[i] = false;
              if (fDebug) std::cout<<"ivtx = "<<i<<" no tracks around the vertex matched to the flash."<<std::endl;
              continue;
          }

          //study if mult>=2
          if (trkindex[i].size()>1){
            //find two longest tracks
            int j0 = -1;
            int j1 = -1;
            double trklen0 = -1;
            double trklen1 = -1;
            //find the highest track
            int jhigh = -1;
            double highy= -1000;
            for (size_t j = 0; j<trkindex[i].size(); ++j){//Loop over all tracks around vertex
              if (trklen[trkindex[i][j]]>trklen0){
                j1 = j0;
                trklen1 = trklen0;
                j0 = j;
                trklen0 = trklen[trkindex[i][j]];
              }
              else if (trklen[trkindex[i][j]]>trklen1){
                j1 = j;
                trklen1 = trklen[trkindex[i][j]];
              }

              if (trkstarty[trkindex[i][j]]>highy){
                highy = trkstarty[trkindex[i][j]];
                jhigh = j;
              }
              if (trkendy[trkindex[i][j]]>highy){
                highy = trkendy[trkindex[i][j]];
                jhigh = j;
              }
            }//Loop over all tracks around vertex
            float dcosx0 = trkstartdcosx[trkindex[i][j0]];
            float dcosy0 = trkstartdcosy[trkindex[i][j0]];
            float dcosz0 = trkstartdcosz[trkindex[i][j0]];
            float dcosx1 = trkstartdcosx[trkindex[i][j1]];
            float dcosy1 = trkstartdcosy[trkindex[i][j1]];
            float dcosz1 = trkstartdcosz[trkindex[i][j1]];
            if (fliptrack[i][j0]){
              dcosx0 = trkenddcosx[trkindex[i][j0]];
              dcosy0 = trkenddcosy[trkindex[i][j0]];
              dcosz0 = trkenddcosz[trkindex[i][j0]];
            }
            if (fliptrack[i][j1]){
              dcosx1 = trkenddcosx[trkindex[i][j1]];
              dcosy1 = trkenddcosy[trkindex[i][j1]];
              dcosz1 = trkenddcosz[trkindex[i][j1]];
            }
            //cosine angle between two longest tracks, 1 indicates broken track
            cosangle[i] = std::abs(dcosx0*dcosx1+dcosy0*dcosy1+dcosz0*dcosz1);
            if (j0 == jhigh){
              dcosylong[i] = fliptrack[i][j0]?std::abs(trkstartdcosy[trkindex[i][j0]]):std::abs(trkenddcosy[trkindex[i][j0]]);
              trklen2nd[i] = trklen[trkindex[i][j1]];
            }
            if (fDebug) std::cout<<i<<" "<<trkindex[i].size()<<" "<<j0<<" "<<jhigh<<" "<<cosangle[i]<<" "<<dcosylong[i]<<" "<<trklen2nd[i]<<std::endl;
            if (cosangle[i]>fMaxCosineAngle) {
              nuvtx[i] = false;
              continue;
            }
            if (j0 == jhigh){
              if (dcosylong[i]>fMaxCosy1stTrk&&trklen2nd[i]<fMinTrackLen2ndTrk){
                nuvtx[i] = false;
                continue;

              }
            }
          }//Multi>1

          //If there are more than 2 tracks, accept the vertex 
          if (trkindex[i].size()>2){
            nuvtx[i] = true;
            if (fDebug) std::cout<<"ivtx = "<<i<<" track multiplicity = "<<trkindex[i].size()<<std::endl;
            continue;
          }

          //Single track
          if (trkindex[i].size()==1){
            nuvtx[i] = false; //will update later
            //Only select contained track
            int itrk = trkindex[i][0];
            //Track is fully contained
            if (inFV(trkstartx[itrk],
                     trkstarty[itrk],
                     trkstartz[itrk])&&
                inFV(trkendx[itrk],
                     trkendy[itrk],
                     trkendz[itrk])){
              if (std::abs(trkstartdcosy[itrk])>fMaxCosySingle) {
                nuvtx[i] = false;
                continue;
              }
              //At least 40 cm
              if (trklen[itrk]>fMinTrackLenSingle){
                double TrackLengthYNu = trklen[itrk]*std::abs(trkstartdcosy[itrk]);
                double LongSingleTrackdEdxRatio = -999;
                if (trkstarty[itrk]>trkendy[itrk]){
                  LongSingleTrackdEdxRatio = trkStartdEdx[itrk]/trkEnddEdx[itrk];
                }
                else{
                  LongSingleTrackdEdxRatio = trkEnddEdx[itrk]/trkStartdEdx[itrk];
                }

                //if (fDebug) std::cout<<"LongSingleTrackdEdxRatio = "<<LongSingleTrackdEdxRatio<<" TrackLengthYNu "<<TrackLengthYNu<<" trkStartdEdx "<<trkStartdEdx[itrk]<<" trkEnddEdx "<<trkEnddEdx[itrk]<<std::endl;
                if(LongSingleTrackdEdxRatio>fMindEdxRatioSingle || (TrackLengthYNu<=fMaxTrkLengthySingle && LongSingleTrackdEdxRatio<=fMindEdxRatioSingle)){
                  nuvtx[i] = true;
                }
              }//At least 40 cm
            }//Track is contained
            if (fDebug) std::cout<<"ivtx = "<<i<<" single track "<<nuvtx[i]<<std::endl;
            continue;
            }//Single track

          //Multiplicity = 2
          if (trkindex[i].size()==2){
            bool isMichel = false;
            double trkstartdedx0 = 0;
            double trkenddedx0 = 0;
            double trkendy0 = 0;
            double trklen1 = 0;
            if (trklen[trkindex[i][0]]>trklen[trkindex[i][1]]){//first track is longer
              if (!fliptrack[i][0]){
                trkstartdedx0 = trkStartdEdx[trkindex[i][0]];
                trkenddedx0 = trkEnddEdx[trkindex[i][0]];
                trkendy0 = trkendy[trkindex[i][0]];
                trklen1 = trklen[trkindex[i][1]];
              }
              else{
                trkstartdedx0 = trkEnddEdx[trkindex[i][0]];
                trkenddedx0 = trkStartdEdx[trkindex[i][0]];
                trkendy0 = trkstarty[trkindex[i][0]];
                trklen1 = trklen[trkindex[i][1]];
              }
            }//first track is longer
            else{//second track is longer
              if (!fliptrack[i][1]){
                trkstartdedx0 = trkStartdEdx[trkindex[i][1]];
                trkenddedx0 = trkEnddEdx[trkindex[i][1]];
                trkendy0 = trkendy[trkindex[i][1]];
                trklen1 = trklen[trkindex[i][0]];
              }
              else{
                trkstartdedx0 = trkEnddEdx[trkindex[i][1]];
                trkenddedx0 = trkStartdEdx[trkindex[i][1]];
                trkendy0 = trkstarty[trkindex[i][1]];
                trklen1 = trklen[trkindex[i][0]];
              }
            }//second track is longer
            if (((trkstartdedx0>trkenddedx0&&
                  trkstartdedx0>fMinStartdEdx1stTrk&&trkenddedx0<fMaxEnddEdx1stTrk)||
                   trkendy0>fDistToEdgeY)&&trklen1<fMinTrackLen2ndTrk)
              isMichel = true;

            if (isMichel){
              nuvtx[i] = false;
            }
            else{
              nuvtx[i] = true;
            }
            if (fDebug) std::cout<<"ivtx = "<<i<<" mul = 2 "<<nuvtx[i]<<std::endl;
            continue;
            }//Multiplicity = 2

        }//Loop over all vertices

        //Find the longest track
        int ivtx = -1;
        int itrk = -1;
        double longesttracklength = -1;
        for (size_t i = 0; i<vtxlist->size(); ++i){
          if (!nuvtx[i]) continue;
          for (size_t j = 0; j<trkindex[i].size(); ++j){
            if (fDebug) std::cout<<"ivtx = "<<i<<" trkid = "<<trkindex[i][j]<<" tracklen "<<trklen[trkindex[i][j]]<<" trackflashmatch "<<trackflashmatch[trkindex[i][j]]<<std::endl;
            if (trklen[trkindex[i][j]]>longesttracklength&&
                trackflashmatch[trkindex[i][j]]&&
                trklen[trkindex[i][j]]>fMinTrackLen){
              longesttracklength = trklen[trkindex[i][j]];
              ivtx = i;
              itrk = trkindex[i][j];
            }
          }
        }//Loop over all vertices
        if (ivtx!=-1 && itrk!=-1){
          if (fDebug) std::cout<<ivtx<<" "<<itrk<<std::endl;
          //outputfile[isample]<<run<<" "<<subrun<<" "<<event<<" "<<ivtx<<" "<<trkindex[ivtx][itrk]<<" "<<trkindex[ivtx].size()<<std::endl;
          //util::CreateAssn(*fMyProducerModule, evt, tracklist[itrk], vtxlist[ivtx], *vertexTrackAssociations);
        }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////  End Selection II Cuts 
    ///////  Begin HitRatio Cut     
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////    


    if ( evt_i != 0 ) genNewHitRatioValues();

    auto vtx = vtxlist->at(ivtx); 

    std::vector<double> vtxXYZ = { vtx.X(), vtx.Y(), vtx.Z() };
    auto vtxWT  = fGeomH->Point_3Dto2D(vtxXYZ,2); //plane);
    
    auto vtx_w = vtxWT.w ; // Comes in cm 
    auto vtx_t = vtxWT.t + 800. * fGeomH->TimeToCm() ; 

    std::vector<::cv::Point> contour;

    //float rad = 60. ;
    std::vector<float> x(50,0);
    std::vector<float> y(50,0);

    int hits_in_rad = 0.; 
    int hits_in_rad_g = 0.; 

    // Build circle around radius
    for(int i = 0; i < x.size(); i++){

       x[i] = vtx_w + fRadius * cos(M_PI * 2. * i / x.size());
       y[i] = vtx_t + fRadius * sin(M_PI * 2. * i / x.size());

       ::cv::Point pt(x[i],y[i]) ;
       contour.emplace_back(pt);
       }

     ::cv::Point pt(contour.at(0).x, contour.at(0).y) ;
     contour.emplace_back(pt);

     float ratio = 0.;

     for(auto const & h : *ev_hit_g){

         if( h.WireID().Plane != 2 ) continue;

         ::cv::Point h_pt(h.WireID().Wire * fGeomH->WireToCm(),
                          h.PeakTime() * fGeomH->TimeToCm()) ;
         auto inside = ::cv::pointPolygonTest(contour,h_pt,false);

         if(inside  >= 0) hits_in_rad_g ++ ;

         if(inside  >= 0 && h.GoodnessOfFit() >= 0) hits_in_rad++ ;
         }

     if( hits_in_rad_g == 0 )
        ratio = 0;
     else
        ratio = float(hits_in_rad)/hits_in_rad_g;

    if( ratio < fHitRatio || hits_in_rad_g < fMinHitsRequired ) return false ;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////  End HitRatio Cuts 
    ///////  Begin Showerreco Cuts     
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////    

    float min_IP = 1e9;
    int min_it = -1 ;

    std::vector<std::pair<int,int>> CandidatePairs;
    std::vector<int> cand_ids;

    for ( int s1 = 0; s1 < ev_s->size(); s1++ ){

            auto const& shr1 = ev_s->at(s1);

            for ( int s2 = 0; s2 < ev_s->size(); s2++ ){

                if (s2 <= s1) continue;

                auto const& shr2 = ev_s->at(s2);

                geoalgo::Vector_t dflipA(-1.*shr1.Direction()) ; 
                geoalgo::Vector_t dflipB(-1.*shr2.Direction()) ;

                // Make the backwards projection for the showers
                auto bksa = ::geoalgo::HalfLine_t(shr1.ShowerStart(),dflipA);
                auto bksb = ::geoalgo::HalfLine_t(shr2.ShowerStart(),dflipB);

                // Calc the Opening angle of the showers
                double Oangle = acos( shr1.Direction().Dot(shr2.Direction())) ; 

                // Calc the vertex point of the two showers. the true designated backwards project
                geoalgo::Point_t vertex(3);
    
                auto st1 = shr1.ShowerStart();
                auto st2 = shr2.ShowerStart();
                auto dir1 = shr1.Direction();
                auto dir2 = shr2.Direction();
                geoalgo::HalfLine_t shr1_hl(st1.X(),st1.Y(),st1.Z(),dir1.X(),dir1.Y(), dir1.Z() );
                geoalgo::HalfLine_t shr2_hl(st2.X(),st2.Y(),st2.Z(),dir2.X(),dir2.Y(), dir2.Z() );

                _geoAlgo.commonOrigin(shr1_hl, shr2_hl, vertex, true);

                // Calc Diretion of two correlated shower
                geoalgo::Vector_t momentum(3);// need to fill out
                geoalgo::Vector_t mom_vect(shr2.Direction()*shr1.Energy(2) +shr1.Direction()*shr2.Energy(2)) ;
                mom_vect.Normalize();
                momentum = mom_vect * sqrt(pow(shr1.Energy(),2)+pow(shr2.Energy(),2)+2*shr2.Energy()*shr1.Energy()*cos(Oangle));

                //===========================================
                auto _RIP = pow(_geoAlgo.SqDist(bksa,bksb),0.5);
                auto _RadL_A = vertex.Dist(shr1.ShowerStart());
                auto _RadL_B = vertex.Dist(shr2.ShowerStart());

                if(Oangle < fMinOpeningAngle ){
                  //std::cout<<"Bad Angle: "<< Oangle<<std::endl ; 
                  continue;
                 }

                if(pow(_geoAlgo.SqDist(bksa,bksb),0.5) > fMaxIP ){
                   //std::cout<<"Ip? "<<pow(_geoAlgo.SqDist(bksa,bksb),0.5)<<std::endl ; 
                   continue;
                  }

                if(_RadL_A > fMaxRadLength || _RadL_B > fMaxRadLength ){
                   //std::cout<<"Rad Length : "<<_RadL_A<<", "<<_RadL_B<<std::endl ;
                   continue;
                    }

                // Bunch of cuts
                if( _RIP < min_IP ){
                  min_IP = _RIP ;
                  min_it = CandidatePairs.size();
                  }

                CandidatePairs.push_back(std::make_pair(s1,s2));
                cand_ids.emplace_back(s1);
                cand_ids.emplace_back(s2);

            }// shower ID 2 
        }// shower ID 1 

        if( CandidatePairs.size() != 1 ) return false; 

        if( isSignal ) _tagged_v[evt_i]++; 
        else _bkgd_v[evt_i]++;

    }// loop over iterations per event


    return true;
  }

bool RepoSelectionII::inFV(double x, double y, double z) const
{
    double distInX = x - fGeometry->DetHalfWidth();
    double distInY = y;
    double distInZ = z - 0.5 * fGeometry->DetLength();

    if (std::abs(distInX) < fDistToEdgeX && std::abs(distInY) < fDistToEdgeY && std::abs(distInZ) < fDistToEdgeZ) return true;

    return false;
}

double RepoSelectionII::scaledEdx(double x, int plane, bool isdata) const{
  double dEdx = 1.63;
  double p0_data[3] = {1.927, 2.215, 1.936};
  double p1_data[3] = {0.001495, 0.0001655, 0.001169};
  double p0_mc[3] = {1.843, 1.904, 1.918};
  double p1_mc[3] = {-0.0008329, -0.001357, -0.0007563};
  if (isdata){
    return dEdx/(p0_data[plane]+x*p1_data[plane]);
  }
  else{
    return  dEdx/(p0_mc[plane]+x*p1_mc[plane]);
  }
}

 double RepoSelectionII::random(double low, double high) {
    double random = ((double) rand()) / (double) RAND_MAX;
    double diff = high - low ;
    return random * diff + low ;
 }

  bool RepoSelectionII::finalize() {

    float FV = 216.35 * 193 * 1016.8 ;  // Fiducial volume 
    float rho = 1.4 ; // g / cm3
    float avogadro = 6.022*pow(10,23);
    float g_per_mole = 39.948 ;
    float n_nucleon = 40 ;
    float n_targ = rho * FV / g_per_mole * avogadro * n_nucleon ;
    float flux = _tot_pot / 1e20 * 1.20e11 ; //1.2e11 is integrated flux from beam group for 1e20 POT , 0.5 - 2.0 GeV

    //float eff = 1. ; // This will be replaced by selection efficiency of loosened Sel2

    for( int i = 0; i < _n_it_per_event; i++ ){
      //eff *= float( _tagged_v[i] - _bkgd_v[i] ) / _n_true_signal ;
      _xsec_v[i] = float(_tagged_v[i] - _bkgd_v[i]) / flux / n_targ ; // / eff ;
      _eff_v[i] = float(_tagged_v[i] - _bkgd_v[i]) / _n_true_signal ;

      }
  
    return true;
  }

}
#endif
