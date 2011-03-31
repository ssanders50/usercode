Double_t centmid[]={2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,45,55,65,75};
Double_t negval[]={0.513,0.685,0.763,0.796,0.804,0.801,0.784,0.755,0.695,0.572,0.405,0.237};
Double_t negerr[]={0.00359,0.00321,0.00318,0.00317,0.00322,0.00319,0.00316,0.00318,0.00224,0.00237,0.00285,0.00424};
Double_t posval[]={0.484, 0.66,0.746,0.782,0.796, 0.79,0.775,0.746,0.689,0.562,0.403,0.247};
Double_t poserr[]={0.0034,0.00309,0.00311,0.00311,0.00318,0.00315,0.00312,0.00314,0.00222,0.00233,0.00283,0.00442};

void RescorPlot(){
  TGraphErrors * Pos = new TGraphErrors(12,centmid,posval,0,poserr);
  TGraphErrors * Neg = new TGraphErrors(12,centmid,negval,0,negerr);
  TH1F * hPos = new TH1F("hPos","",500,0,100);
  hPos->SetMaximum(1.0);
  hPos->SetStats(kFALSE);
  Pos->SetMarkerStyle(20);
  Pos->SetMarkerSize(1.2);
  Pos->SetMarkerColor(kBlue);
  Pos->SetLineColor(kBlue);
  hPos->SetXTitle("Centrality");
  hPos->SetYTitle("EP Resolution Correction Factor");
  Neg->SetMarkerStyle(21);
  Neg->SetMarkerSize(1.2);
  Neg->SetMarkerColor(kRed);
  Neg->SetLineColor(kRed);
  Neg->SetLineColor(kRed);
  TLegend * leg = new TLegend(0.6,0.68,0.86,0.88);
  leg->AddEntry(Neg,"-2 #leq #eta < -1","p");
  leg->AddEntry(Pos,"  1 #leq #eta < 2","p");
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  hPos->GetYaxis()->SetTitleOffset(0.9);
  hPos->Draw();
  Pos->Draw("p");
  Neg->Draw("p");
  leg->Draw("p");
  c1->Print("RescorPlot.pdf","pdf");
  c1->Print("RescorPlot.png","png");

}
