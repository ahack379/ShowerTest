//Author: Anne Schukraft

int AS_Flux() {

   //double targetpot = 5.283e+19 * 20000/19750;
   double targetpot = 2.42e20;

   //BNB flux
   TFile *flux = new TFile("numode_bnb_470m_r200.root");
   TH1D *h_flux_numu = (TH1D*) flux -> Get("numu");

   //TFile *flux = new TFile("./fluxfiles/numode_bnb_470m_r200_from_simple_ntuples.root");
   //TH1D *h_flux_numu = (TH1D*) flux -> Get("hsum");
   double fluxpot = 1e20;

   h_flux_numu -> Scale(targetpot/fluxpot);
   cout << "The total flux for " << targetpot << " POT is " << h_flux_numu -> Integral() << endl;

   TCanvas *c = new TCanvas("c", "c", 800, 600);
   c -> SetTicks(1, 1);
   c -> SetGridx(1);
   c -> SetGridy(1);
   h_flux_numu -> SetStats(0);
   (h_flux_numu -> GetYaxis()) -> SetTitleOffset(1.5);
   h_flux_numu -> Draw();

   TH1F *h_flux_numu_cut = (TH1F*)h_flux_numu->Clone("h_flux_numu_cut");
   for(int i = 0; i < h_flux_numu_cut -> GetNbinsX(); i++) {
      if(h_flux_numu_cut -> GetBinCenter(i) < 0.){
        std::cout<<"SHIT "<<std::endl ;
        h_flux_numu_cut -> SetBinContent(i, 0);
	}
      //if(h_flux_numu_cut -> GetBinCenter(i) < 0.4 || h_flux_numu_cut->GetBinCenter(i) > 4.) h_flux_numu_cut -> SetBinContent(i, 0);
   }
   //h_flux_numu_cut -> SetMarkerStyle(20);
   //h_flux_numu_cut -> SetMarkerColor(4);
   h_flux_numu_cut -> SetLineWidth(4);
   h_flux_numu_cut -> Draw("same");
  
   double mean = h_flux_numu_cut-> GetMean();
   cout << "The mean energy is: " << mean << endl;
   int binmean = h_flux_numu_cut -> FindBin(mean);
   cout << "The bin of the mean is: " << binmean << endl;

   int n = h_flux_numu_cut -> GetNbinsX();

   double lowerint = h_flux_numu_cut -> Integral(1, binmean);
   cout << lowerint << endl;
   double lowerborder = lowerint * 0.32;
   cout << lowerborder << endl;
   double lowersum = 0;
   int i = 0;
   while (lowersum < lowerborder) {
      i++;
      lowersum += h_flux_numu_cut -> GetBinContent(i);
      cout << i << "\t" << lowersum << endl;
   }

   cout << "The total flux for " << targetpot << " POT is " << h_flux_numu -> Integral() << endl;
   cout << "The total flux for " << targetpot << " POT is " << h_flux_numu_cut -> Integral() << endl;

   cout << lowersum << endl;
   double low = h_flux_numu_cut -> GetBinCenter(i-1);
   cout << "The lower edge bin is: " << i-1 << endl;
   cout << "The lower edge center energy is: " << low << endl;
   cout << "The lower energy error is: " << mean - low << endl;

   double upperint = h_flux_numu_cut -> Integral(binmean, n);
   cout << upperint << endl;
   double upperborder = upperint * 0.32;
   double uppersum = 0;
   i = 0;
   while (uppersum < upperborder) {
      uppersum += h_flux_numu_cut -> GetBinContent(n+1 - i);
      i++;
   }

   double up = h_flux_numu_cut -> GetBinCenter(n+1 - (i-1));
   cout << "The upper edge bin is: " << i-1 << endl;
   cout << "The upper edge center energy is: " << up << endl;
   cout << "The upper energy error is: " << up - mean << endl;

   TGraph *gmean = new TGraph();
   gmean -> SetPoint(0, mean, 0);
   gmean -> SetPoint(1, mean, 1e10);
   gmean -> SetLineWidth(2);
   gmean -> SetLineColor(kOrange+1);
   gmean -> Draw("same");

   //TGraph *glow = new TGraph();
   //glow -> SetPoint(0, low, 0);
   //glow -> SetPoint(1, low, 1e10);
   //glow -> SetLineWidth(2);
   //glow -> SetLineColor(kOrange+1);
   //glow -> SetLineStyle(7);
   //glow -> Draw("same");

   //TGraph *gup = new TGraph();
   //gup -> SetPoint(0, up, 0);
   //gup -> SetPoint(1, up, 1e10);
   //gup -> SetLineWidth(2);
   //gup -> SetLineColor(kOrange+1);
   //gup -> SetLineStyle(7);
   //gup -> Draw("same");

   TLegend *l = new TLegend(0.6, 0.7, 0.89, 0.89);
   l -> AddEntry(h_flux_numu, "BNB #nu_{#mu} flux, #nu-mode", "l");
   l -> AddEntry(gmean, "<E_{#nu}>", "l");
   //l -> AddEntry(glow, "1#sigma energy range", "l");
   l -> Draw();

   return 0;
}
