//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class larlite::CCpi0Eff+;
#pragma link C++ class larlite::CrossSection+;
#pragma link C++ class larlite::BackgroundCalc+;
#pragma link C++ class larlite::SingleGammaBackgroundCalc+;
//#pragma link C++ class larlite::POTCalc+;
#pragma link C++ class larlite::AcceptanceStudy+;
#pragma link C++ class larlite::BackgroundAll+;
#pragma link C++ class larlite::BackgroundObnoxious+;
#pragma link C++ class larlite::CCEff+;
#pragma link C++ class larlite::NuMode+;
#pragma link C++ class larlite::MCNuClusterer+;
#pragma link C++ class larlite::BackgroundBT+;
#pragma link C++ class larlite::AngleBackgroundBT+;
#pragma link C++ class larlite::PostSel2ShowerQuality+;
//#pragma link C++ class larlite::MCTrackClusterer+;

//ADD_NEW_CLASS ... do not change this line
#endif

