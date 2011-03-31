//#include "HiEvtPlaneList.h"
#include "../../../RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include <iomanip>
#include <fstream>
using namespace std;
TFile * tf;

void makeRescor(){

  TString tag = "_hiGoodMergedTracks_03_18PtCut_PtWeight";
  // Step 1.   Generate EP Resolutions
  tf = new TFile(Form("../data/rpflat_combined%s.root",tag.Data()));
  TH1D * c12[50];
  TH1D * c13[50];
  TH1D * c23[50];
  TH1D * ccnt12[50];
  TH1D * ccnt13[50];
  TH1D * ccnt23[50];
  TH1D * rescor[50];
  
  for(int i = 0; i< NumEPNames; i++) {
    TString baseName = Form("v2analyzer/v2/v2Reco/%s",EPNames[i].data());
    c12[i] = (TH1D *) tf->Get(Form("%s/c12_%s",baseName.Data(),EPNames[i].data()));
    c13[i] = (TH1D *) tf->Get(Form("%s/c13_%s",baseName.Data(),EPNames[i].data()));
    c23[i] = (TH1D *) tf->Get(Form("%s/c23_%s",baseName.Data(),EPNames[i].data()));
    ccnt12[i] = (TH1D *) tf->Get(Form("%s/ccnt12_%s",baseName.Data(),EPNames[i].data()));
    ccnt13[i] = (TH1D *) tf->Get(Form("%s/ccnt13_%s",baseName.Data(),EPNames[i].data()));
    ccnt23[i] = (TH1D *) tf->Get(Form("%s/ccnt23_%s",baseName.Data(),EPNames[i].data()));
    c12[i]->Divide(ccnt12[i]);
    c13[i]->Divide(ccnt13[i]);
    c23[i]->Divide(ccnt23[i]);
    rescor[i] = (TH1D *) c12[i]->Clone(Form("rescor_%s",EPNames[i].data()));
    rescor[i]->Multiply(c13[i]);
    rescor[i]->Divide(c23[i]);
    rescor[i]->SetTitle(EPNames[i].data());
    for(int j = 1; j<=c12[i]->GetNbinsX(); j++) {
      if(rescor[i]->GetBinContent(j)>0) {
	rescor[i]->SetBinContent(j, TMath::Sqrt(rescor[i]->GetBinContent(j)));
      } else {
	rescor[i]->SetBinContent(j,0);
      }
    }
  }
  TFile * fRescor = new TFile(Form("rescor/%s.root",tag.Data()),"RECREATE");
  fRescor->cd();
  for(int i = 0; i<NumEPNames; i++) {
    rescor[i]->Write();
  }
  fRescor->Close();
}

