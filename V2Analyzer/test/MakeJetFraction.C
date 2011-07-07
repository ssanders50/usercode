TH2D * cntAll = 0;
TH2D * cntJet = 0;
TH2D * ratio = 0;
void MakeJetFraction() {
  TFile * tf = new TFile("../data/rpflat_combined.root");
  TH2D * cnt0 = 0;
  TH2D * cnt1 = 0;
  TH2D * cnt2 = 0;
  TH2D * cnt3 = 0;
  TH2D * cnt4 = 0;
  for(ieta = 21; ieta<= 30; ieta++) {
    if(ieta==21) {
      cntAll = (TH2D *) tf->Get("v2analyzer/v2/v2Reco/etHF/cnt_etHF_21");
      cnt0 =   (TH2D *) tf->Get("v2analyzer/v2/v2Reco/etHF/jet0cnt_etHF_21");
      cnt1 =   (TH2D *) tf->Get("v2analyzer/v2/v2Reco/etHF/jet1cnt_etHF_21");
      cnt2 =   (TH2D *) tf->Get("v2analyzer/v2/v2Reco/etHF/jet2cnt_etHF_21");
      cnt3 =   (TH2D *) tf->Get("v2analyzer/v2/v2Reco/etHF/jet3cnt_etHF_21");
      cnt4 =   (TH2D *) tf->Get("v2analyzer/v2/v2Reco/etHF/jet4cnt_etHF_21");
    } else {
      cntAll->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/etHF/cnt_etHF_%d",ieta)));
      cnt0->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/etHF/jet0cnt_etHF_%d",ieta)));
      cnt1->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/etHF/jet1cnt_etHF_%d",ieta)));
      cnt2->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/etHF/jet2cnt_etHF_%d",ieta)));
      cnt3->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/etHF/jet3cnt_etHF_%d",ieta)));
      cnt4->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/etHF/jet4cnt_etHF_%d",ieta)));
    }
  }
  cntJet = (TH2D *) cnt2->Clone("cntJet");
  cntJet->Add(cnt3);
  cntJet->Add(cnt4);
  ratio = (TH2D *) cntJet->Clone("ratio");
  ratio->Divide(cntAll);
  TCanvas * cratio = new TCanvas("TracksWithJetRatio","TracksWithJetRatio",800,600);
  ratio->SetStats(kFALSE);
  cratio->SetRightMargin(0.2);
  ratio->Draw();
  cratio->Print("TracksWithJetRatio.png","png");
}
