TFile * tf;
void FlattenExample(){
  tf = new TFile("../data/rpflat_combined.root");
  TH1D * hRaw = (TH1D *) tf->Get("hiEvtPlaneFlat/EvtPTracksNegEtaGap/psi");
  TH1D * hFlat = (TH1D *) tf->Get("hiEvtPlaneFlat/EvtPTracksNegEtaGap/psiFlat");
  hRaw->SetStats(kFALSE);
  hFlat->SetLineColor(kRed);
  hRaw->SetTitle("");
  hRaw->SetMaximum(10000);
  hRaw->SetXTitle("#Psi_{2}");
  hRaw->SetYTitle("Counts");
  hRaw->Draw();
  hFlat->Draw("same");
  TLatex * lab = new TLatex(-2,9000,"Tracking Event Plane (-2 #leq #eta < -1)");
  lab->Draw();
  c1->Print("FlattenExample.pdf","pdf");
}
