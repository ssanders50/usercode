void NpartBinCompare(){
  TFile * tf95 = new TFile("../data/rpflat_combined_hiGoodMergedTracks_Flow_Skim_Run2010-v5_eff95.root");
  TFile * tf97 = new TFile("../data/rpflat_combined_hiGoodMergedTracks_Flow_Skim_Run2010-v5.root");
  TFile * tf100 = new TFile("../data/rpflat_combined_hiGoodMergedTracks_Flow_Skim_Run2010-v5_eff100.root");
  TH1D * NpartBin95 = tf95->Get("v2analyzer/NpartBin");
  TH1D * NpartBin97 = tf97->Get("v2analyzer/NpartBin");
  TH1D * NpartBin100 = tf100->Get("v2analyzer/NpartBin");
  TH1D * NpartBinCnt95 = tf95->Get("v2analyzer/NpartBinCnt");
  TH1D * NpartBinCnt97 = tf97->Get("v2analyzer/NpartBinCnt");
  TH1D * NpartBinCnt100 = tf100->Get("v2analyzer/NpartBinCnt");
  NpartBin95->Divide(NpartBinCnt95);
  NpartBin97->Divide(NpartBinCnt97);
  NpartBin100->Divide(NpartBinCnt100);
  NpartBin95->SetMarkerColor(kRed);
  NpartBin100->SetMarkerColor(kBlue);
  NpartBin95->Divide(NpartBin97);
  NpartBin100->Divide(NpartBin97);
  NpartBin95->SetMinimum(0.7);
  NpartBin95->SetMaximum(1.3);
  NpartBin95->Draw();
  NpartBin100->Draw("same");
  TLegend * leg = new TLegend(0.2,0.2,0.5,0.4);
  leg->AddEntry(NpartBin95,"err95/default","lp");
  leg->AddEntry(NpartBin100,"err100/default","lp");
  leg->Draw();
}
