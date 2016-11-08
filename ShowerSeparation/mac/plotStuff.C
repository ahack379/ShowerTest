void plotStuff(){

   TFile * f = new TFile("ana.root","READ");
   TTree * t = (TTree *)f->Get("tree");
   
   std::cout<<"Entries "<< t->GetEntries();
   
   TCanvas * c1 = new TCanvas("c1","c1");
   TH1D * h1 = new TH1D("h1","h1",20,100000,1000000);
   TH1D * h2 = new TH1D("h2","h2",20,100000,1000000);
   t->Draw("sum_charge00>>h1","");
   t->Draw("sum_charge02>>h2","");
   h1->SetFillColorAlpha(kRed,0.35); 
   h2->SetFillColorAlpha(kBlue,0.35); 
   //h4->SetFillStyle( 3001);
   //h5->SetFillStyle( 3001);
   //h5->SetLineColor(kRed);
   h1->GetXaxis()->SetTitle("Charge Integral [ADC]");
   h2->DrawCopy() ;
   h1->Draw("same");

}
