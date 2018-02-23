#ifndef LARLITE_SAVEWEIGHTS_CXX
#define LARLITE_SAVEWEIGHTS_CXX

#include "SaveWeights.h"
#include "DataFormat/vertex.h"
#include "DataFormat/mceventweight.h"
#include <sstream>

namespace larlite {

  bool SaveWeights::initialize() {

    _func_v = {"AGKYpT_Genie","AGKYxF_Genie","DISAth_Genie","DISBth_Genie","DISCv1u_Genie","DISCv2u_Genie","FermiGasModelKf_Genie", "FermiGasModelSf_Genie","FormZone_Genie", "IntraNukeNabs_Genie", "IntraNukeNcex_Genie", "IntraNukeNel_Genie", "IntraNukeNinel_Genie", "IntraNukeNmfp_Genie", "IntraNukeNpi_Genie", "IntraNukePIabs_Genie", "IntraNukePIcex_Genie", "IntraNukePIel_Genie", "IntraNukePIinel_Genie", "IntraNukePImfp_Genie", "IntraNukePIpi_Genie", "NC_Genie", "NonResRvbarp1pi_Genie", "NonResRvbarppi_Genie", "NonResRvp1pi_Genie", "NonResRvppi_Genie", "ResDecayEta_Genie", "ResDecayGamma_Genie", "ResDecayTheta_Genie", "ccresAxial_Genie", "ccresVector_Genie", "cohMA_Genie", "cohR0_Genie", "ncelAxial_Genie", "ncelEta_Genie", "ncresAxial_Genie", "ncresVector_Genie", "qema_Genie", "qevec_Genie"}; 

    if ( _event_producer == "fluxeventweight" )
      _func_v = {"SkinEffect", "HornCurrent", "K-", "K+", "K0", "NucleonInXsec", "NucleonQEXsec", "NucleonTotXsec", "pi-", "piInelasticXsec",  "piQEXsec",  "piTotalXsec",  "pi+"} ;

    _file.open("weights_all.txt",std::ios_base::in);
    std::vector<double> wgt_v ; 
    int i = 0;

    for(std::string line; std::getline(_file, line); )   //read stream line by line
    {   
        std::istringstream in(line);      //make a stream for the line itself
    
        float run;
        in >> run;                  //and read the first whitespace-separated token
        float subrun;
        in >> subrun;                  //and read the first whitespace-separated token
        float event ;
        in >> event ;

        std::vector<double> wgt_v ; 

        auto evinfo = std::make_pair(subrun,event);
        double wgt ;
        while( in >> wgt ){
          wgt_v.emplace_back(wgt);
        }
        //std::cout<<"LENGTH :"<<wgt_v.size()<<", "<<i<<std::endl ;
	    i++;
        _wgtmap[evinfo] = wgt_v ;
    }   

    std::cout<<"Function size: " <<_func_v.size()<<std::endl;

    return true;
  }
  
  bool SaveWeights::analyze(storage_manager* storage) {

    auto run = storage->run_id();
    auto event = storage->event_id();
    auto subrun = storage->subrun_id();

    //auto ev_wgt = storage->get_data<event_mceventweight>("genieeventweight");
    auto ev_wgt = storage->get_data<event_mceventweight>(_event_producer);

    storage->set_id(run, subrun, event);

    std::pair<int,int> evinfo;
    evinfo = std::make_pair(subrun,event);

    // loop through all entries in map
    for (auto const& element : _wgtmap) {

	if ( element.first != evinfo) continue;
	
	//std::cout<<"Comparing subrun + event : "<<element.first.first <<", "<<element.first.second<<" with "<<evinfo.first<<", "<<evinfo.second<<std::endl;

	auto const& wgt_v = element.second;
	std::map<std::string,std::vector<double>> w_map; 

    //if ( _event_producer == "genieeventweight"){

    //      for( int i = 0; i < _func_v.size() ; i++ ){
    //        std::vector<double> temp_wgt_v = {wgt_v[2*i], wgt_v[2*i + 1] };
    //        //std::cout<<"FUNCTION: "<<_func_v[i]<<", "<<temp_wgt_v[0]<<", "<<temp_wgt_v[1]<<std::endl ;

    //        w_map[_func_v[i]] = temp_wgt_v ; 
    //      }

    //   larlite::mceventweight thisweight(w_map);
    //   ev_wgt->emplace_back( thisweight);

    //   if (ev_wgt->size() == 0 )
    //     std::cout<<"WHAT IS HAPPENING\n\n\n\n\n\n\n\n "<<std::endl ;

    //      return true;
    //    
    // }
    // else{

	  for( int i = 0; i < _func_v.size() ; i++ ){
	    std::vector<double> temp_wgt_v; 
        for ( int j = 0; j < 1000; j++ )
          temp_wgt_v.emplace_back(wgt_v[1000*i + j]);

	    w_map[_func_v[i]] = temp_wgt_v ; 
	  }

       larlite::mceventweight thisweight(w_map);
       ev_wgt->emplace_back( thisweight);

       if (ev_wgt->size() == 0 )
         std::cout<<"WHAT IS HAPPENING\n\n\n\n\n\n\n\n "<<std::endl ;

	  return true;
	
    //  }
    }

    return true;
  }

  bool SaveWeights::finalize() {

    return true;
  }

}
#endif
