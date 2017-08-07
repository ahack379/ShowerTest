//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class larlite::AVFilter+;
#pragma link C++ class larlite::FlashCut+;
#pragma link C++ class larlite::RatioCut+;
#pragma link C++ class larlite::EnergyFilter+;


#pragma link C++ class larlite::TrackMultiplicity+;
#pragma link C++ class larlite::Pi0Cuts+;
#pragma link C++ class larlite::SingleShowerPi0Cuts+;
#pragma link C++ class larlite::DistanceCut+;
//ADD_NEW_CLASS ... do not change this line
#endif


