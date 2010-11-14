#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TList.h"
#include "TString.h"
#include "TFile.h"
#include "TPaveText.h"
#include <time.h>
#include <cstdlib>

#include "../../../RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"

using namespace std;
#include <vector>
using std::vector;
using std::rand;

TFile * tf;
void
checkResolution()
{
  using namespace std;
  tf = new TFile("../data/rpflat_combined.root");
  //
  //Get Event Planes
  //
  TH1D * c12[50];
  TH1D * c13[50];
  TH1D * c23[50];
  TH1D * ccnt12[50];
  TH1D * ccnt13[50];
  TH1D * ccnt23[50];
  TH1D * rescor[50];
  TString tag = "Hydjet_Bass";
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
  TCanvas * can[8];
  TLegend * leg[8];
  TPaveText * basedon[50];

  Int_t indx=0;
  for(int i = 0; i<NumEPNames; i=i+8) {
    TString cname = Form("EP_%d_%d",i,i+7);
    can[indx] = new TCanvas(cname.Data(),cname.Data(),1150,800);
    can[indx]->Divide(4,2);
    for(int j = 0; j< 8; j++) {
      if(i+j>=NumEPNames) break;
      can[indx]->cd(j+1);
      rescor[i+j]->SetMaximum(1.2);
      rescor[i+j]->SetXTitle("centrality");
      rescor[i+j]->SetYTitle("ResCor");
      rescor[i+j]->Draw(); 
      rescor[i+j]->SetStats(kFALSE);
      can[indx]->Update();
      TPaveText * tp = (TPaveText *) gPad->GetPrimitive("title");
      tp->SetX2NDC(1.8*tp->GetX2NDC());
      tp->SetY1NDC(tp->GetY2NDC() - 1.8*(tp->GetY2NDC()-tp->GetY1NDC()));
      tp->Draw();
      tp->SetTextSize(0.06);
      basedon[i+j] = new TPaveText(20.,0.9,95.,1.15);
      TString name1 = EPNames[RCMate1[i+j]];
      TString name2 = EPNames[RCMate2[i+j]];
      cout<<name1.Data()<<" "<<name2.Data()<<endl;
      basedon[i+j]->AddText(Form("%5.1f < #eta < %5.1f",EPEtaMin1[i+j],EPEtaMax1[i+j]));
      if(EPEtaMin2[i+j] != 0 || EPEtaMax2[i+j]!=0) {
	basedon[i+j]->AddText(Form("%5.1f < #eta < %5.1f",EPEtaMin2[i+j],EPEtaMax2[i+j]));
      }
      basedon[i+j]->AddText("Based on: ");
      basedon[i+j]->AddText(Form("%s (%5.1f < #eta < %5.1f)",name1.Data(),EPEtaMin1[RCMate1[i+j]],EPEtaMax1[RCMate1[i+j]]));
      basedon[i+j]->AddText(Form("%s (%5.1f < #eta < %5.1f)",name2.Data(),EPEtaMin1[RCMate2[i+j]],EPEtaMax1[RCMate2[i+j]]));
      basedon[i+j]->Draw();
    }
    can[indx]->Print(Form("~/public_html/ResCor/%s_%s.png",tag.Data(),can[indx]->GetName()),"png");
    ++indx;
  }
}

