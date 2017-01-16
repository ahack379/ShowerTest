#include <string>
//Author: Anne Schukraft

using std::string;
using namespace std;

int CCpizero() {

   //energy (not momentum!) thresholds in GeV
   double Ar_Ethresh_p = 0;
   double Ar_Ethresh_mu = 0;
   double Ar_Ethresh_pi = 0;
   double Ar_Ethresh_pi0 = 0;
   double Ar_Ethresh_K = 0;
   double Ar_Ethresh_K0 = 0;
   double Ar_Ethresh_e = 0;
   double Ar_Ethresh_Lambda0 = 0;
   double Ar_Ethresh_Sigma0 = 0;
   double Ar_Ethresh_Sigma = 0;
   double Ar_Ethresh_antip = 0;

   double C_Ethresh_p = 0;
   double C_Ethresh_mu = 0;
   double C_Ethresh_pi = 0;
   double C_Ethresh_pi0 = 0;
   double C_Ethresh_K = 0;
   double C_Ethresh_K0 = 0;
   double C_Ethresh_e = 0;
   double C_Ethresh_Lambda0 = 0;
   double C_Ethresh_Sigma0 = 0;
   double C_Ethresh_Sigma = 0;
   double C_Ethresh_antip = 0;

   double nargon = 1;
   double ncarbon = 1;

   int Ar_n = 17;
   double Ar_E[] = {100, 150, 200, 250, 300, 400, 600, 800, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 10000, 30000};
   double *Ar_sigma = new double[Ar_n];

   TFile *Ar_f_numu = new TFile("/Users/ah673/WorkArea/Root6LArLite/UserDev/FindNeutrinos/CalcEfficiency/mac/GenieFiles/gxsec_numu_argon.root");
   TDirectory *Ar_dir_numu = (TDirectory*) Ar_f_numu -> Get("nu_mu_Ar40");
   TGraph *Ar_genie_cc_numu = (TGraph*) Ar_dir_numu -> Get("tot_cc");

   //std::cout<<"GRAPH STUFF: "<<Ar_genie_cc_numu->GetN();
   for(int i = 0; i < Ar_n; i++) {
      Ar_sigma[i] = Ar_genie_cc_numu -> Eval(Ar_E[i]/1000.);
    //  std::cout<<Ar_sigma[i]<<", "<<Ar_E[i]/1000<<std::endl;
   }

   int C_n = 10;
   double C_E[] = {150, 200, 250, 300, 400, 600, 800, 1000, 1500, 2000};
   double *C_sigma = new double[C_n];

   TFile *C_f_numu = new TFile("/Users/ah673/WorkArea/Root6LArLite/UserDev/FindNeutrinos/CalcEfficiency/mac/GenieFiles/gxsec_numu_carbon.root");
   TDirectory *C_dir_numu = (TDirectory*) C_f_numu -> Get("nu_mu_C12");
   TGraph *C_genie_cc_numu = (TGraph*) C_dir_numu -> Get("tot_cc");
   for(int i = 0; i < C_n; i++) {
      C_sigma[i] = C_genie_cc_numu -> Eval(C_E[i]/1000.);
   }

   double *Ar_N_ccincl = new double[Ar_n];
   double *Ar_N_0pi = new double[Ar_n];
   double *C_N_ccincl = new double[C_n];
   double *C_N_0pi = new double[C_n];

   int nf;
   int *pdgf = new int[1000];
   double *Ef = new double[1000]; 

   //loop numu - argon files
   char filename[200];
   cout << "Reading numu - argon!" << endl;

   for(int k = 0; k < Ar_n; k++) {

      double wght = 0;

      cout << "... " << k << endl;

      int energy = Ar_E[k];
      TTree *tree = 0;
      sprintf(filename, "/Users/ah673/WorkArea/Root6LArLite/UserDev/FindNeutrinos/CalcEfficiency/mac/GenieFiles/gntpc_numu_argon_%dMeV_CCincl.root", energy);
      TFile file(filename,"READ");

      tree = dynamic_cast <TTree *>(file.Get("gst"));
      if(!tree) return 1;

      tree -> SetBranchAddress("nf", &nf);
      tree -> SetBranchAddress("pdgf", pdgf);
      tree -> SetBranchAddress("Ef", Ef);

      int nev = tree -> GetEntries();
      //std::cout<<"Nev: "<<nev<<std::endl ;

      double n_ccincl = 0;
      double n_0pi = 0;

      int npi = 0;
      int npi0 = 0;
      int np = 0;
      int nk = 0;
      int nk0 = 0;
      int ngamma = 0;
      int nn = 0;
      int npseudo = 0;
      int nlambda0 = 0;
      int nsigma0 = 0;
      int nsigma = 0;
      int nantip = 0;
      int nexclude = 0;

      for(int i = 0; i < nev; i++) {

         tree->GetEntry(i);

         npi = 0;
         npi0 = 0;
         np = 0;
         nk = 0;
         nk0 = 0;
         ngamma = 0;
         nn = 0;
         nlambda0 = 0;
         nsigma0 = 0;
         nsigma = 0;
         nantip = 0;
         nexclude = -1;
         nelec = 0;
         ////std::cout<<"\n";
         for(int j = 0; j < nf; j++) {
	    //std::cout<<"PDG: "<<pdgf[j]<<", ";
            if(pdgf[j] == 111 && Ef[j] > Ar_Ethresh_pi0) npi0++;
            else if(pdgf[j] == 11 )  nelec++;
            else if((pdgf[j] == 211 || pdgf[j] == -211) && Ef[j] > Ar_Ethresh_pi) npi++;
            else if(pdgf[j] == 2212 && Ef[j] > Ar_Ethresh_p) np++;
            else if((pdgf[j] == 321 || pdgf[j] == -321) && Ef[k] > Ar_Ethresh_K) nk++;
            else if((pdgf[j] == 130 || pdgf[j] == 310 || pdgf[j] == 311 || pdgf[j] == -311) && Ef[j] > Ar_Ethresh_K0) nk0++;
            else if(pdgf[j] == 22) ngamma++;
            else if(pdgf[j] == 2112 || pdgf[j] == -2112) nn++;
            else if((pdgf[j] == 3122 || pdgf[j] == -3122) && Ef[j] > Ar_Ethresh_Lambda0) nlambda0++;
            else if(pdgf[j] == 3212 && Ef[j] > Ar_Ethresh_Sigma0) nsigma0++;
            else if((pdgf[j] == 3222 || pdgf[j] == 3112) && Ef[j] > Ar_Ethresh_Sigma) nsigma++;
            else if(pdgf[j] == -2212 && Ef[j] > Ar_Ethresh_antip) nantip++;
         }// end loop over particles	

         n_ccincl++;
         if(npi0 == 1 && npi == 0 && nk == 0 && nk0 == 0 && ngamma == 0 && nelec == 0) n_0pi++;

      }//end loop over events

      file.Close();
   
      Ar_N_ccincl[k] = n_ccincl / nev / nargon * Ar_sigma[k];
      Ar_N_0pi[k] = n_0pi / nev / nargon * Ar_sigma[k];
   }


   for(int k = 0; k < C_n; k++) {

      double wght = 0;

      cout << "... " << k << endl;

      int energy = C_E[k];
      TTree *tree = 0;
      sprintf(filename, "/Users/ah673/WorkArea/Root6LArLite/UserDev/FindNeutrinos/CalcEfficiency/mac/GenieFiles/gntpc_numu_carbon_%dMeV_CCincl.root", energy);
      TFile file(filename,"READ");

      tree = dynamic_cast <TTree *>(file.Get("gst"));
      if(!tree) return 1;

      tree -> SetBranchAddress("nf", &nf);
      tree -> SetBranchAddress("pdgf", pdgf);
      tree -> SetBranchAddress("Ef", Ef);

      int nev = tree -> GetEntries();

      double n_ccincl = 0;
      double n_0p = 0;
      double n_1p = 0;
      double n_2p = 0;
      double n_0pi = 0;

      int npi = 0;
      int npi0 = 0;
      int np = 0;
      int nk = 0;
      int nk0 = 0;
      int ngamma = 0;
      int nn = 0;
      int npseudo = 0;
      int nlambda0 = 0;
      int nsigma0 = 0;
      int nsigma = 0;
      int nantip = 0;
      int nexclude = 0;

      for(int i = 0; i < nev; i++) {

         tree->GetEntry(i);

         npi = 0;
         npi0 = 0;
         np = 0;
         nk = 0;
         nk0 = 0;
         ngamma = 0;
         nelec = 0;
         nn = 0;
         nlambda0 = 0;
         nsigma0 = 0;
         nsigma = 0;
         nantip = 0;
         nexclude = -1;

         for(int j = 0; j < nf; j++) {
            if(pdgf[j] == 111 && Ef[j] > C_Ethresh_pi0) npi0++;
	    else if(abs(pdgf[j]) == 11 ) nelec++;
            else if((pdgf[j] == 211 || pdgf[j] == -211) && Ef[j] > C_Ethresh_pi) npi++;
            else if(pdgf[j] == 2212 && Ef[j] > C_Ethresh_p) np++;
            else if((pdgf[j] == 321 || pdgf[j] == -321) && Ef[k] > C_Ethresh_K) nk++;
            else if((pdgf[j] == 130 || pdgf[j] == 310 || pdgf[j] == 311 || pdgf[j] == -311) && Ef[j] > C_Ethresh_K0) nk0++;
            else if(pdgf[j] == 22) ngamma++;
            else if(pdgf[j] == 2112 || pdgf[j] == -2112) nn++;
            else if((pdgf[j] == 3122 || pdgf[j] == -3122) && Ef[j] > C_Ethresh_Lambda0) nlambda0++;
            else if(pdgf[j] == 3212 && Ef[j] > C_Ethresh_Sigma0) nsigma0++;
            else if((pdgf[j] == 3222 || pdgf[j] == 3112) && Ef[j] > C_Ethresh_Sigma) nsigma++;
            else if(pdgf[j] == -2212 && Ef[j] > C_Ethresh_antip) nantip++;
         }// end loop over particles    


         nexclude = npi + npi0 + nk + nk0 + ngamma + nlambda0 + nsigma + nsigma0 + nantip;
         n_ccincl++;
         //if(npi0 == 1) n_0pi++;
         if(npi0 == 1 && npi == 0 && nk == 0 && nk0 == 0 && ngamma == 0 && nelec == 0) n_0pi++;
         //if(npi0 == 1 && (npi != 0 || nk != 0 || nk0 != 0 || ngamma != 0 || nelec != 0) ) 
	   //std::cout<<"WELL...OK..."<<std::endl ;

      }//end loop over events

      file.Close();

      C_N_ccincl[k] = n_ccincl / nev / ncarbon * C_sigma[k];
      C_N_0pi[k] = n_0pi / nev / ncarbon * C_sigma[k];
   }




   TCanvas *c = new TCanvas("c", "c", 1000, 600);
   c->SetGrid();
   //c -> SetLogx(1);

   TGraph *fake = new TGraph();
   fake -> SetPoint(0, 100, 0);
   //fake -> SetPoint(1, 2000,3); //100000, 100);
   fake -> SetPoint(1, 3000,3); //100000, 100);
   fake -> SetMarkerColor(10);
   (fake -> GetXaxis()) -> SetTitle("E_{#nu} [MeV]");
   (fake -> GetYaxis()) -> SetTitle("#sigma (10^{-38} cm^2 / nucleon");

   fake -> Draw("AP");
   TGraph *gr_Ar_0pi = new TGraph(Ar_n, Ar_E, Ar_N_0pi);
   gr_Ar_0pi -> SetTitle("#nu - Ar");
   gr_Ar_0pi -> SetLineWidth(2);
   gr_Ar_0pi -> SetLineStyle(1);
   (gr_Ar_0pi -> GetXaxis()) -> SetRangeUser(10, 2000);//100000);
   (gr_Ar_0pi -> GetYaxis()) -> SetRangeUser(0, 2);
   (gr_Ar_0pi -> GetXaxis()) -> SetTitle("E_{#nu} [MeV]");
   (gr_Ar_0pi -> GetYaxis()) -> SetTitle("#sigma (10^{-38} cm^2 / nucleon");

   auto x = gr_Ar_0pi->GetX() ;
   auto y = gr_Ar_0pi->GetY() ;
   float area = 0;
   for ( int i = 0; i < 16; i++){
     area += (y[i] * (x[i+1] - x[i]) ) ; // /2000);
     std::cout<<"X,Y: "<< x[i]<<", "<<y[i]<<std::endl; 
     }


     std::cout<<"Total area and ave: "<<area<<", "<<area/(2000 - 100)<<std::endl ;

   //std::cout<<"Integral! "<<gr_Ar_0pi->Integral(0,-1)<<", "<<gr_Ar_0pi->Integral()/2000 <<std::endl ;

   gr_Ar_0pi -> Draw("LSAME");
   TGraph *gr_C_0pi = new TGraph(C_n, C_E, C_N_0pi);
   gr_C_0pi -> SetLineWidth(2);
   gr_C_0pi -> SetLineStyle(1);
   gr_C_0pi -> SetLineColor(2);
   gr_C_0pi -> Draw("LSAME");
   //TGraph *gr_Ar_incl = new TGraph(Ar_n, Ar_E, Ar_N_ccincl);
   //gr_Ar_incl -> SetTitle("#nu - Ar");
   //gr_Ar_incl -> SetLineWidth(2);
   //gr_Ar_incl -> SetLineStyle(7);
   //gr_Ar_incl -> Draw("LSAME");
   //TGraph *gr_C_incl = new TGraph(C_n, C_E, C_N_ccincl);
   //gr_C_incl -> SetLineWidth(2);
   //gr_C_incl -> SetLineStyle(7);
   //gr_C_incl -> SetLineColor(2);
   //gr_C_incl -> Draw("LSAME");

   cout << "Done!" << endl;

   return 0;
}

