
#include "/net/hisrv0001/home/sanders/CMSSW_3_9_2/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include <iomanip>
#include <fstream>
using namespace std;
TFile * tf;

void checkRotations(){

  tf = new TFile("$CMSSW_BASE/src/V2Analyzer/V2Analyzer/data/rpflat_combined.root");
  TH2D * h2d[5];
  TH1D * h1low[5];
  TH1D * h1high[5];
  h2d[0] = (TH2D *) tf->Get("v2analyzer/TRACKp_etHFp");
  h2d[1] = (TH2D *) tf->Get("v2analyzer/TRACKp_HFp");
  h2d[4] = (TH2D *) tf->Get("v2analyzer/TRACKm_HFm");
  h2d[3] = (TH2D *) tf->Get("v2analyzer/TRACKm_etHFm");
  h2d[2] = (TH2D *) tf->Get("v2analyzer/TRACKp_TRACKm");
  TCanvas * c1 = new TCanvas("c1","c1",700,550);
  c1->Divide(3,2);
  for(int i = 0; i<5; i++) {
    c1->cd(i+1);
    gPad->SetGrid(1,1);
    h2d[i]->SetStats(kFALSE);
    h2d[i]->Draw();
    h1low[i] = (TH1D *) h2d[i]->ProjectionX(Form("low%d",i),15,20);
    h1high[i] = (TH1D *) h2d[i]->ProjectionX(Form("high%d",i),80,85);
    h1high[i]->SetLineColor(2);
    h1low[i]->SetStats(kFALSE);
    h1high[i]->SetStats(kFALSE);
  }
  c1->Print("~/public_html/ReactPlaneCorrelations2D.png");
  TCanvas * c2 = new TCanvas("c2","c2",700,550);
  c2->Divide(3,2);
  for(int i = 0; i<5; i++) {
    c2->cd(i+1);
    gPad->SetGrid(1,1);
    h1low[i]->Draw();
    h1high[i]->Draw("same");
  }
  c2->Print("~/public_html/ReactPlaneCorrelations1D.png");

}

