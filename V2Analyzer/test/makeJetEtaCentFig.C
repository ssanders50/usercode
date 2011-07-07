
void makeJetEtaCentFig() {
  TFile * tf = new TFile("../data/rpflat_combined.root");
  TH2D * jetEtaCent = (TH2D *) tf ->Get("v2analyzer/jetEtaCent");
  TH2D * nojetEtaCent = (TH2D *) tf->Get("v2analyzer/nojetEtaCent");
  jetEtaCent->SetStats(kFALSE);
  jetEtaCent->SetOption("colz");
  nojetEtaCent->SetStats(kFALSE);
  nojetEtaCent->SetOption("colz");
  TCanvas * c = new TCanvas("JetEtaCent","JetEtaCent",850,600);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetRightMargin(0.15);
  jetEtaCent->SetXTitle("#eta (Jet p_{T}#geq30)");
  jetEtaCent->SetYTitle("Centrality");
  jetEtaCent->Draw();
  c->cd(2);
  gPad->SetRightMargin(0.15);
  nojetEtaCent->Draw();
  nojetEtaCent->SetXTitle("#eta (Jet p_{t}<30)");
  nojetEtaCent->SetYTitle("Centrality");
  c->Print("JetEtaCent.png","png");
}
