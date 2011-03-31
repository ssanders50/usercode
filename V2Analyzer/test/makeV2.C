//#include "HiEvtPlaneList.h"
#include "../../../RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include <iomanip>
#include <fstream>
using namespace std;
TFile * tf;
TFile * fout;
static const Int_t nCentBins = 14;
static const Int_t nPtBins = 15;
//static const Int_t nEtBins = 15;
//static const Int_t nEtaBins = 22;
static const Int_t nEtaBins = 50;
static const double minpt = 0.3;
static const double maxpt = 3.0;
static const double centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,90,100};
static const double ptbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0};
static const double etbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0};
//static const double etabins[]={-5,-4.5,-4,-3.5,-3,-2.4,-2,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,3,3.5,4,4.5,5};
static const double etabins[]={-5.0, -4.8, -4.6, -4.4, -4.2, -4.0, -3.8, -3.6, -3.4, -3.2,
			       -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2,
			       -1.0, -0.8, -0.6, -0.4, -0.2,  0.0,  0.2,  0.4,  0.6,  0.8,
                               1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,  2.8,
			        3.0,  3.2,  3.4,  3.6,  3.8,  4.0,  4.2,  4.4,  4.6,  4.8,
			        5.0};
// STANDARD ETA BINS:     |eta|<0.8
//                  0.8 < |eta| < 1.6
//                  1.6 < |eta| < 2.0
//                  2.0 < |eta| < 2.4
TH1D * hcentbins = new TH1D("centbins","centbins",nCentBins,centbins);
TH1D * hptbins = new TH1D("ptbins","ptbins",nPtBins,ptbins);
TH1D * hetabins = new TH1D("hetabins","hetabins",nEtaBins,etabins);

void AddToV2(Int_t type, Int_t rp, int minc, int maxc,  Double_t mineta, Double_t maxeta, TH1D ** rescor);
double IntegralV2(Int_t icent,TH1D * pt,TGraphErrors * g,TGraphErrors * gpt, double &err, double & dndeta, double & dndetaerr);
TGraphErrors * GenV2(double centmin, double centmax, TH2D * hpt, Int_t marker, Int_t color);
void makeV2Prog(Int_t EP1, double mineta1, double maxeta1, 
	    Int_t EP2, double mineta2,   double maxeta2,
	    TString tag);
TH2D * EffCorr(TH2D * hPt, TH2D * hPtCnt, Double_t eta);

TH2D * EffCorr(TH2D * hPt, TH2D * hPtCnt, Double_t eta){
  TString effName = Form("eff_%05.1",eta);
  TH2D * eff = (TH2D *) hPtCnt->Clone(effName.Data());
  for(int j = 1; j<=eff->GetNbinsY(); j++) {
    Double_t cent = eff->GetYaxis()->GetBinCenter(j); 
    Int_t ictab = (cent+0.1)/5;
    TH2D * heffRef = (TH2D *) feff->Get(Form("rEff_cbin%d",ictab));
    TH2D * hfakRef = (TH2D *) feff->Get(Form("rFak_cbin%d",ictab));
    for(int i = 1; i<=eff->GetNbinsX(); i++) {
      Double_t pt = hPt->GetBinContent(i,j);
      Double_t ptcnt = hPtCnt->GetBinContent(i,j);
      if(ptcnt>0) {
	pt/=ptcnt;
      } else {
	continue;
      }
      Int_t iptab = heffRef->GetYaxis()->FindBin(pt);
      Int_t iptmin = iptab;
      Int_t iptmax = iptab;
      Double_t ptmin = 0;
      Double_t ptmax = 0;
      if(pt>=heffRef->GetYaxis()->GetBinCenter(iptab)) {
	iptmax = iptmin+1;
	ptmin = heffRef->GetYaxis()->GetBinCenter(iptmin);
	ptmax = heffRef->GetYaxis()->GetBinCenter(iptmax);
      } else {
	iptmin=iptmax-1;
	ptmin = heffRef->GetYaxis()->GetBinCenter(iptmin);
	ptmax = heffRef->GetYaxis()->GetBinCenter(iptmax);
      }

      Int_t ietatab = heffRef->GetXaxis()->FindBin(eta);
      Double_t ef = heffRef->GetBinContent(ietatab,iptab);
      Double_t efmin = heffRef->GetBinContent(ietatab,iptmin);
      Double_t efmax = heffRef->GetBinContent(ietatab,iptmax);
      ef = (pt-ptmin)*(efmax-efmin)/(ptmax-ptmin)+efmin;
      Double_t fa = hfakRef->GetBinContent(ietatab,iptab);
      Double_t famin = hfakRef->GetBinContent(ietatab,iptmin);
      Double_t famax = hfakRef->GetBinContent(ietatab,iptmax);
      fa = (pt-ptmin)*(famax-famin)/(ptmax-ptmin)+famin;
      if(fa>0.99) fa = 0.99;
      eff->SetBinContent(i,j,ef/(1-fa));
      eff->SetBinError(i,j,0.);
    }
  }
  return eff;
}

Bool_t GenCalc = kFALSE;
TH2D * cos1_;
TH2D * cnt1_;
TH2D * cos_;
TH2D * cnt_;
TFile * feff = 0;
Double_t MINETA1;
Double_t MAXETA1;
Double_t MINETA2;
Double_t MAXETA2;
Double_t MINCENT;
Double_t MAXCENT;
void makeV2() {
  //TString tag = "Flow_Skim_Run2010-v5_dz5Flat_-10to10";
  TString tag = "dzerr14_chi280";

  if(minpt!=0.3) tag=tag+Form("-minpt%d",(Int_t) (10.*minpt+0.001));
  if(maxpt!=3.0) tag=tag+Form("-maxpt%d",(Int_t) (10.*maxpt+0.001));
  cout<<"tag: "<<tag.Data()<<endl;
  fout = new TFile(Form("EPSpectra_%s.root",tag.Data()),"RECREATE");
  fout->cd();
  TH1D * ptbins = new TH1D("ptbins","ptbins",nPtBins,ptbins);
  fout->Write();
  makeV2Prog(EvtPTracksPosEtaGap, -0.8,  0.0, EvtPTracksNegEtaGap, 0.0, 0.8, tag);
  //makeV2Prog(EvtPTracksPosEtaGap, -2.4, -2.0, EvtPTracksNegEtaGap, 2.0, 2.4, tag);
  //makeV2Prog(EvtPTracksPosEtaGap, -2.0, -1.6, EvtPTracksNegEtaGap, 1.6, 2.0, tag);
  //makeV2Prog(EvtPTracksPosEtaGap, -1.6, -1.2, EvtPTracksNegEtaGap, 1.2, 1.6, tag);
  //makeV2Prog(EvtPTracksPosEtaGap, -1.2, -0.8, EvtPTracksNegEtaGap, 0.8, 1.2, tag);
  //makeV2Prog(EvtPTracksPosEtaGap, -0.8, -0.4, EvtPTracksNegEtaGap, 0.4, 0.8, tag);
  //makeV2Prog(EvtPTracksPosEtaGap, -0.4,  0.0, EvtPTracksNegEtaGap, 0.0, 0.4, tag);
  makeV2Prog(EvtPTracksPosEtaGap, -2.4, -2.0, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksPosEtaGap, -2.0, -1.6, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksPosEtaGap, -1.6, -1.2, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksPosEtaGap, -1.2, -0.8, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksPosEtaGap, -0.8, -0.4, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksPosEtaGap, -0.4,  0.0, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksNegEtaGap,  0.0,  0.4, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksNegEtaGap,  0.4,  0.8, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksNegEtaGap,  0.8,  1.2, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksNegEtaGap,  1.2,  1.6, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksNegEtaGap,  1.6,  2.0, -1, 0, 0, tag);
  makeV2Prog(EvtPTracksNegEtaGap,  2.0,  2.4, -1, 0, 0, tag);
}
void makeV2Prog(Int_t EP1 , double mineta1 , double maxeta1 , 
		Int_t EP2 , double mineta2 ,   double maxeta2 ,
		TString tag){

  //--------------------------------------
  // Set up analysis
  //
  Int_t type = 1;  // =1 for tracks, =2 for etcalo, = 3 for hcal(special replay)
  double aveta1 = (fabs(mineta1)+fabs(maxeta1))/2.;
  cout<<EPNames[EP1].data()<<" "<<mineta1<<" "<<maxeta1<<endl;
  if(EP2>=0) cout<<EPNames[EP2].data()<<" "<<mineta2<<" "<<maxeta2<<endl;
  TString effFileName ="trkCorrFlow_hydjet100k_badWedge_SQ_5cent_vertexZ10.root";
  //TString effFileName="trkCorrFlow_hydjet100k_badWedge_SQ_5cent_dz10chi40_vertexZ10.root";
  if(!effFileName.Contains("NoEffCorr")) feff = new TFile(effFileName.Data());
  if(GenCalc) tag = tag+"_GENERATOR";
  Int_t MaxCentBin = 12;
  //
  // End of setup
  //-------------------------------------
  MINETA1 = mineta1;
  MAXETA1 = maxeta1;
  MINETA2 = mineta2;
  MAXETA2 = maxeta2;
  Double_t deta1 = maxeta1-mineta1;
  Double_t deta2 = maxeta2-mineta2;
  cos_ = 0;
  cnt_ = 0;
  cos1_ = 0;
  cnt1_ = 0;
  Int_t ietamin1 =  hetabins->FindBin(mineta1);
  Int_t ietamax1 = hetabins->FindBin(maxeta1-0.1);
  Int_t ietamin2 = 0;
  Int_t ietamax2 = 0;
  TString EtaBins = Form("etabins_%d_%d",ietamin1,ietamax1);
  mineta1 = hetabins->GetBinLowEdge(ietamin1);
  maxeta1 = hetabins->GetBinLowEdge(ietamax1)+hetabins->GetBinWidth(ietamax1);
  if(mineta2<maxeta2) {
    ietamin2 = hetabins->FindBin(mineta2);
    ietamax2 = hetabins->FindBin(maxeta2-0.1);
    EtaBins += Form("_%d_%d",ietamin2,ietamax2);
    mineta2 = hetabins->GetBinLowEdge(ietamin2);
    maxeta2 = hetabins->GetBinLowEdge(ietamax2)+hetabins->GetBinWidth(ietamax2);
  }
  Int_t nEta1 = 0;
  Int_t inceta1[20];
  Int_t nEta2 = 0;
  Int_t inceta2[20];
  for(int i = ietamin1; i<=ietamax1; i++) inceta1[nEta1++] = i;
  for(int i = 0; i< nEta1; i++) cout<<"Include eta bin (1): "<<inceta1[i]<<endl;
  if(ietamin2>0) {
    for(int i = ietamin2; i<=ietamax2; i++) inceta2[nEta2++] = i;
    for(int i = 0; i< nEta2; i++) cout<<"Include eta bin (2): "<<inceta2[i]<<endl;
  } 
  cout<<"Eta ranges: "<<mineta1<<"-"<<maxeta1<<" :  "<<mineta2<<"-"<<maxeta2<<endl;
  //                   0-5   5-10  10-15    15-20    20-25    25-30     30-35   35-40    40-50  50-60   60-70     70-80 80-90 90-100
  int mincent[14] ={     0,     5,    10,      15,      20,      25,      30,     35,       40,    50,     60,       70,   80,    90};
  int maxcent[14] ={     5,    10,    15,      20,      25,      30,      35,     40,       50,    60,     70,       80,   90,   100};
  int markers[14]= {    22,    21,    23,      24,      25,      30,      20,      3,       29,    28,     26,       27,    3,     5};
  int colors[14] = {kBlack, kBlue,  kRed, kCyan+2, kViolet, kSpring, kOrange, kRed+4, kAzure+9, kTeal, kRed-4,kYellow+2,    0,     0};

  tf = new TFile("../data/rpflat_combined.root");
  TH1D * hcent = (TH1D *) tf->Get("v2analyzer/cent");

  // Step 0.  Generate Npart is spectra

  TH2D * hpt1[20];
  TH2D * dNdPt1[20];
  TH2D * het1[20];
  TH2D * dNdEt1[20];
  for(int i = 0; i< nEta1; i++) {
    hpt1[i] = (TH2D *) tf->Get(Form("v2analyzer/Spectra/pt_%d",inceta1[i]-1));
    dNdPt1[i] = (TH2D *) tf->Get(Form("v2analyzer/Spectra/ptCnt_%d",inceta1[i]-1));
    het1[i] = (TH2D *) tf->Get(Form("v2analyzer/Spectra/et_%d",inceta1[i]-1));
    dNdEt1[i] = (TH2D *) tf->Get(Form("v2analyzer/Spectra/etCnt_%d",inceta1[i]-1));
    TH2D * efc = EffCorr(hpt1[i],dNdPt1[i], hetabins->GetBinCenter(inceta1[i]));
    hpt1[i]->Divide(efc);
    dNdPt1[i]->Divide(efc);
  }
  
  TH2D * hpt2[20];
  TH2D * dNdPt2[20];
  TH2D * het2[20];
  TH2D * dNdEt2[20];
  if(ietamin2>0) {
    for(int i = 0; i< nEta2; i++) {
      hpt2[i] = (TH2D *) tf->Get(Form("v2analyzer/Spectra/pt_%d",inceta2[i]-1));
      dNdPt2[i] = (TH2D *) tf->Get(Form("v2analyzer/Spectra/ptCnt_%d",inceta2[i]-1));
      het2[i] = (TH2D *) tf->Get(Form("v2analyzer/Spectra/et_%d",inceta2[i]-1));
      dNdEt2[i] = (TH2D *) tf->Get(Form("v2analyzer/Spectra/etCnt_%d",inceta2[i]-1));
      TH2D * efc = EffCorr(hpt2[i],dNdPt2[i], hetabins->GetBinCenter(inceta2[i]));
      hpt2[i]->Divide(efc);
      dNdPt2[i]->Divide(efc);
      
    }
  }

  TH2D * dNdPt1S = (TH2D *) dNdPt1[0]->Clone("dNdPt1S");
  TH2D * dNdEt1S = (TH2D *) dNdEt1[0]->Clone("dNdEt1S");
  TH2D * hpt1S = (TH2D *) hpt1[0]->Clone("hpt1S");
  TH2D * het1S = (TH2D *) het1[0]->Clone("het1S");

  for(int i = 1; i<nEta1; i++) {
    dNdPt1S->Add(dNdPt1[i]);
    dNdEt1S->Add(dNdEt1[i]);
    hpt1S->Add(hpt1[i]);
    het1S->Add(het1[i]);
  }
  hpt1S->Divide(dNdPt1S);
  het1S->Divide(dNdEt1S);
  if(deta1!=deta2) cout<<"WARNING:  ETA Assymmetry! with deta/deta2= "<<deta1<<"/"<<deta2<<endl;
  for(int i = 1; i<= dNdPt1S->GetNbinsX(); i++ ) {
    for(int j = 1; j<=dNdPt1S->GetNbinsY(); j++) {
      dNdPt1S->SetBinContent(i,j, dNdPt1S->GetBinContent(i,j)/
			    (dNdPt1S->GetXaxis()->GetBinWidth(i)*deta1));
    }
  }
  for(int i = 1; i<= dNdEt1S->GetNbinsX(); i++ ) {
    for(int j = 1; j<=dNdEt1S->GetNbinsY(); j++) {
      dNdEt1S->SetBinContent(i,j, dNdEt1S->GetBinContent(i,j)/
			    (dNdEt1S->GetXaxis()->GetBinWidth(i)*deta1));
    }
  }
  TH2D * dNdPt2S = 0;
  TH2D * dNdEt2S = 0;
  TH2D * hpt2S = 0;
  TH2D * het2S = 0;
  if(ietamin2>0) {
    dNdPt2S = (TH2D *) dNdPt2[0]->Clone("dNdPt2S");
    dNdEt2S = (TH2D *) dNdEt2[0]->Clone("dNdEt2S");
    hpt2S = (TH2D *) hpt2[0]->Clone("hpt2S");
    het2S = (TH2D *) het2[0]->Clone("het2S");
    for(int i = 1; i<nEta2; i++) {
      dNdPt2S->Add(dNdPt2[i]);
      dNdEt2S->Add(dNdEt2[i]);
      hpt2S->Add(hpt2[i]);
      het2S->Add(het2[i]);
    }
    hpt2S->Divide(dNdPt2S);
    het2S->Divide(dNdEt2S);
    for(int i = 1; i<= dNdPt2S->GetNbinsX(); i++ ) {
      for(int j = 1; j<=dNdPt2S->GetNbinsY(); j++) {
	dNdPt2S->SetBinContent(i,j, dNdPt2S->GetBinContent(i,j)/
			       (dNdPt2S->GetXaxis()->GetBinWidth(i)*deta2));
      }
    }
    for(int i = 1; i<= dNdEt2S->GetNbinsX(); i++ ) {
      for(int j = 1; j<=dNdEt2S->GetNbinsY(); j++) {
	dNdEt2S->SetBinContent(i,j, dNdEt2S->GetBinContent(i,j)/
			       (dNdEt2S->GetXaxis()->GetBinWidth(i)*deta2));
      }
    }
  }
  TH2D * hpt = (TH2D *) hpt1S->Clone("hpt");
  TH2D * het = (TH2D *) het1S->Clone("het");
  TH2D * dNdPt = (TH2D *) dNdPt1S->Clone("dNdPt");
  TH2D * dNdEt = (TH2D *) dNdEt1S->Clone("dNdEt");
  if(ietamin2>0) {
    hpt->Add(hpt2S);
    het->Add(het2S);
    dNdPt->Add(dNdPt2S);
    dNdEt->Add(dNdEt2S);
    hpt->Scale(0.5);
    het->Scale(0.5);
    dNdPt->Scale(0.5);
    dNdEt->Scale(0.5);
  }
  TH1D * NpartBin = (TH1D *) tf->Get("v2analyzer/NpartBin");
  TH1D * NpartBinCnt = (TH1D *) tf->Get("v2analyzer/NpartBinCnt");
  NpartBin->Divide(NpartBinCnt);
  // Step 1.   Generate EP Resolutions
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
    //rescor[i]->Sumw2();
    rescor[i]->Multiply(c13[i]);
    rescor[i]->Divide(c23[i]);
    rescor[i]->SetTitle(EPNames[i].data());
    for(int j = 1; j<=c12[i]->GetNbinsX(); j++) {
      if(rescor[i]->GetBinContent(j)>0) {
	double val = rescor[i]->GetBinContent(j);
	double err = rescor[i]->GetBinError(j);
	rescor[i]->SetBinContent(j, TMath::Sqrt(val));
	rescor[i]->SetBinError(j,0.5*err/TMath::Sqrt(val));
      } else {
	rescor[i]->SetBinContent(j,0);
	rescor[i]->SetBinError(j,0);
      }
    }
  }
  TFile * fRescor = new TFile(Form("rescor/%s.root",tag.Data()),"RECREATE");
  fRescor->cd();
  for(int i = 0; i<NumEPNames; i++) {
    rescor[i]->Write();
    //cout<<EPNames[i].data()<<endl;
    //for(int j = 0; j<rescor[i]->GetNbinsX(); j++) {
    //cout<<mincent[j]<<"\t"<<maxcent[j]<<"\t"<<rescor[i]->GetBinContent(j+1)<<"\t"<<rescor[i]->GetBinError(j+1)<<endl;
    //}
  }
  fRescor->Close();
  tf->cd();
  TCanvas * cRescor;
  if(mineta2==maxeta2) {
    cRescor = new TCanvas(Form("Rescor_%s",EPNames[EP1].data()),Form("Rescor_%s",EPNames[EP1].data()),400,500);
    rescor[EP1]->Draw();
  } else {
    Bool_t special = kTRUE;
    if(!special) {
      cRescor = new TCanvas(Form("Rescor_%s_%s",EPNames[EP1].data(),EPNames[EP2].data()),Form("Rescor_%s_%s",EPNames[EP1].data(),EPNames[EP2].data()),650,550);
    } else {
      cRescor = new TCanvas(Form("Rescor_%s_%s_special",EPNames[EP1].data(),EPNames[EP2].data()),Form("Rescor_%s_%s",EPNames[EP1].data(),EPNames[EP2].data()),850,550);
      cRescor->Divide(2);
      cRescor->cd(1);
    }
    rescor[EP1]->SetMarkerStyle(22);
    rescor[EP1]->SetMarkerColor(kBlue);
    rescor[EP1]->SetMarkerSize(1.5);
    rescor[EP2]->SetMarkerStyle(23);
    rescor[EP2]->SetMarkerColor(kMagenta);
    rescor[EP2]->SetMarkerSize(1.5);
    rescor[EP1]->SetXTitle("Centrality");
    rescor[EP1]->SetYTitle("R (Resolution Correction Factor)");
    gPad->SetGrid(1,1);
    rescor[EP1]->SetMaximum(1.0);
    rescor[EP1]->SetMinimum(0.);
    rescor[EP1]->Draw();
    rescor[EP2]->Draw("same");
    TLegend * legRescor = new TLegend(0.40,0.8,0.92,0.92);
    legRescor->AddEntry(rescor[EP1],EPNames[EP1].data(),"p");
    legRescor->AddEntry(rescor[EP2],EPNames[EP2].data(),"p");
    legRescor->Draw();
    if(special) {
      cRescor->cd(2);
      rescor[EPTracksPosEtaBigGap]->SetMarkerStyle(22);
      rescor[EPTracksPosEtaBigGap]->SetMarkerColor(kBlue);
      rescor[EPTracksPosEtaBigGap]->SetMarkerSize(1.5);
      rescor[EPTracksNegEtaBigGap]->SetMarkerStyle(23);
      rescor[EPTracksNegEtaBigGap]->SetMarkerColor(kMagenta);
      rescor[EPTracksNegEtaBigGap]->SetMarkerSize(1.5);
      rescor[EPTracksPosEtaBigGap]->SetXTitle("Centrality");
      rescor[EPTracksPosEtaBigGap]->SetYTitle("R (Resolution Correction Factor)");
      gPad->SetGrid(1,1);
      rescor[EPTracksPosEtaBigGap]->SetMaximum(1.0);
      rescor[EPTracksPosEtaBigGap]->SetMinimum(0.);
      rescor[EPTracksPosEtaBigGap]->Draw();
      rescor[EPTracksNegEtaBigGap]->Draw("same");
    TLegend * legRescor2 = new TLegend(0.40,0.8,0.92,0.92);
      legRescor2->AddEntry(rescor[EPTracksPosEtaBigGap],EPNames[EPTracksPosEtaBigGap].data(),"p");
      legRescor2->AddEntry(rescor[EPTracksNegEtaBigGap],EPNames[EPTracksNegEtaBigGap].data(),"p");
      legRescor2->Draw();
    }
    cRescor->Print(Form("rescor/%s_%s.pdf",cRescor->GetName(),tag.Data()),"pdf");

  }
  Double_t xx[25];
  Double_t yy[25];
  Double_t xxerr[25];
  Double_t yyerr[25];
  TCanvas * c2 = new TCanvas("c2","c2",800,600);
  gPad->SetGrid(1,1);

  TH1D * hfig;
  if(type == 1 ) {
    hfig = new TH1D("hfig",Form("%s",tag.Data()),100,0,12);
  } else if (type == 2 || type == 3) {
    hfig = new TH1D("hfig",Form("%s",tag.Data()),100,0,6);
  }
  hfig->SetMinimum(0.0);
  hfig->SetMaximum(0.3);
  if(type==2) {
    hfig->SetMinimum(-0.4);
    hfig->SetMaximum(0.6);
  }
  hfig->SetStats(kFALSE);
  TString typeName;
  if(type==1) {
    hfig->SetXTitle("p_{T} (GeV/c)");
    typeName = "Track";
  } else if (type==2) {
    hfig->SetXTitle("E_{T} (GeV)");
    typeName = "Calo";
  } else if (type == 3) {
    hfig->SetXTitle("E_{T} (GeV)");
    typeName = "HCal";
  }
  hfig->SetYTitle("v_{2}");
  hfig->Draw();
  TPaveText * desc;
  if(type==1) {
    desc = new TPaveText(0.1,0.25,3.6,0.29);
  } else if (type == 2 || type==3) {
    desc = new TPaveText(0.1,0.25,3.8,0.29);
  }
  TString label = Form("RP: %s",EPNames[EP1]);
  if(GenCalc) label += "_GENERATOR";
  desc->AddText(label.Data());

  if(type==1) {
    desc->AddText(Form("%5.1f #leq #eta_{track} < %5.1f",mineta1,maxeta1));
  } else if (type==2) {
      desc->AddText(Form("%5.1f #leq #eta_{calo} < %5.1f",mineta1,maxeta1));
  } else if (type==3) {
      desc->AddText(Form("%5.1f #leq #eta_{hcal} < %5.1f",mineta1,maxeta1));
  }
  if(mineta2<maxeta2) {
    desc->AddText(Form("RP: %s",EPNames[EP2]));
    if(type==1) {
      desc->AddText(Form("%5.1f #leq #eta_{track} < %5.1f",mineta2,maxeta2));
    } else if (type==2) {
      desc->AddText(Form("%5.1f #leq #eta_{calo} < %5.1f",mineta2,maxeta2));
    } else if (type==3) {
      desc->AddText(Form("%5.1f #leq #eta_{hcal} < %5.1f",mineta2,maxeta2));
    }
  }
  desc->Draw();

  TH1D * pt[12];
  TGraphErrors * gpt[12];
  TLegend * leg2 = new TLegend(0.60,0.50,0.85,0.88,"Centrality    (N_{part})");
  for(int icent = 0; icent<MaxCentBin; icent++) {
    if(type==1) {
      pt[icent] = (TH1D *) dNdPt->ProjectionX(Form("hpt_%d_%d",mincent[icent],maxcent[icent]),icent+1,icent+1);
    } else { 
      pt[icent] = (TH1D *) dNdEt->ProjectionX(Form("het_%d_%d",mincent[icent],maxcent[icent]),icent+1,icent+1);
    }
    int nptbins = pt[icent]->GetNbinsX();
    int ncnt = 0;
    for(int i = 0; i<nptbins; i++ ) {
      if(pt[icent]->GetBinContent(i+1)>0) {
	Double_t scale = 1./hcent->Integral(hcent->FindBin(mincent[icent]),hcent->FindBin(maxcent[icent]-0.1));
	yy[ncnt]= scale*pt[icent]->GetBinContent(i+1);
	yyerr[ncnt]=scale*pt[icent]->GetBinError(i+1);
	if(yy[ncnt]>0&&pt[icent]->GetBinError(i+1)<sqrt(pt[icent]->GetBinContent(i+1))) yyerr[ncnt] = scale*sqrt(pt[icent]->GetBinContent(i+1));
	xxerr[ncnt]=0;
	if(type==1) {
	  xx[ncnt]=hpt->GetBinContent(i+1,icent+1);
	} else if (type==2 || type==3) {
	  xx[ncnt]=het->GetBinContent(i+1,icent+1);
	}
	++ncnt;
      } 
      
    }
    gpt[icent] = new TGraphErrors(ncnt,xx,yy,xxerr,yyerr);
    gpt[icent]->SetMarkerStyle(markers[icent]);
    gpt[icent]->SetMarkerColor(colors[icent]);
    gpt[icent]->SetLineColor(colors[icent]);
    //gpt[icent]->Draw("p");
    leg2->AddEntry(gpt[icent],  Form("%d-%d     (%5.1f)",mincent[icent],maxcent[icent],NpartBin->GetBinContent(icent+1)),"lp");
  }

  TGraphErrors * g[12];
  TLegend * leg = new TLegend(0.75,0.67,0.92,0.92,"Centrality    (N_{part})");
  std::ofstream file;
  double emax = fabs(mineta1);
  if(fabs(maxeta1)>emax) emax = fabs(maxeta1);
  for(Int_t i = 0; i<MaxCentBin; i++) {
    TString baseName = Form("EP_%d-%d_%02d",mincent[i],maxcent[i],(Int_t)(10.*emax));
    TString fileName =; 
    if(deta1!=deta2) {
      if(mineta2==maxeta2 && maxeta1<=0) {
	baseName+="_NegEta";
      } else if(mineta2==maxeta2 && mineta1>=0) {
	baseName+="_PosEta";
      } else {
	baseName+="_Asym";
      }
    }
    if(tag.Length()>0) {
      baseName+="-";
      baseName+=tag.Data();
      baseName+="_";
      baseName+=EPNames[EP1];
      if(mineta2!=maxeta2) {
	baseName+="_";
	baseName+=EPNames[EP2];
      }
    }
    TString fileName = "results/v2pt_";
    fileName+=baseName;
    fileName+=".txt";
    file.open(fileName.Data());
    if(cos_) {
      cos_->Reset();
      cnt_->Reset();
    }
    if(cos1_) {
      cos1_->Reset();
      cnt1_->Reset();
    }
    AddToV2(type,EP1,mincent[i],maxcent[i],mineta1,maxeta1,rescor); 
    if(mineta2<maxeta2) {
      cos1_ = (TH2D *) cos_->Clone("cos1_");
      cnt1_ = (TH2D *) cnt_->Clone("cnt1_");
      cos_->Reset();
      cnt_->Reset();
      AddToV2(type,EP2,mincent[i],maxcent[i],mineta2,maxeta2,rescor);
      cos_->Add(cos1_);
      cnt_->Add(cnt1_);
    }
    if(type == 1 ) {
      g[i]   = GenV2(mincent[i],maxcent[i],hpt, markers[i],  colors[i]);
    } else if (type == 2||type==3) {
      g[i]   = GenV2(mincent[i],maxcent[i],het, markers[i],  colors[i]);
    }
    g[i]->Draw("p");
    leg->AddEntry(g[i],  Form( "%d-%d     (%5.1f)",mincent[i],maxcent[i],NpartBin->GetBinContent(i+1)),"lp");
    file<<Form("<pt>\tv2\tv2 err\t(1/Nevt)d2N/dptdeta\terr\t%Centrality(d-%d)     (Npart = %5.1f)  ",mincent[i],maxcent[i],NpartBin->GetBinContent(i+1))<<endl;
    Double_t * xxx = g[i]->GetX();
    Double_t * yyy = g[i]->GetY();  
    Double_t * yyyerr = g[i]->GetEY();
    for(int j = 0; j< g[i]->GetN();  j++ ) {
      if(yyy[j]<0.0001||xxx[j]<0.3) continue;
      double * dndeta = gpt[i]->GetY();
      double * dndetaerr = gpt[i]->GetEY();
      file<<setprecision(3)<<xxx[j]<<"\t"<<setprecision(3)<<yyy[j]<<"\t"<<yyyerr[j]<<"\t"<<setprecision(3)<<dndeta[j]<<"\t"<<setprecision(3)<<dndetaerr[j]<<endl;
    }
    file<<tag.Data()<<endl;
    file<<endl;
    file.close();
  } 
  leg->Draw();
  if(mineta2==maxeta2) {
    c2->Print(Form("~/public_html/V2Spectra/v2_%s_%s_%s_%s.png",typeName.Data(),EPNames[EP1].data(),EtaBins.Data(),tag.Data()),"png");
  } else {
    c2->Print(Form("~/public_html/V2Spectra/v2_%s_%s_%s_%s_%s.png",typeName.Data(),EPNames[EP1].data(),EPNames[EP2].data(),EtaBins.Data(),tag.Data()),"png");
  }
  TCanvas * c3; 
  if(type==1 ) {
    c3= new TCanvas("ptDist","ptDist",800,600);
  } else if (type==2||type==3) {
    c3= new TCanvas("etDist","etDist",800,600);
  }
  TH1D * hptframe;
  if(type == 1 ) {
    hptframe = new TH1D("hptframe",tag.Data(),100,0,9);
  } else if (type == 2 || type==3) {
    hptframe = new TH1D("hptframe",tag.Data(),100,0,12);
  }
  hptframe->SetMaximum(10000);
  hptframe->SetMinimum(0.001);
  hptframe->SetStats(kFALSE);
  gPad->SetLogy();
  if(type==1) {
    hptframe->SetXTitle("p_{T} (GeV/c)"); 
    hptframe->SetYTitle("#frac{1}{N_{evt}}#frac{d^{2}N}{dp_{T}d#eta }");
  } else if(type==2||type==3) {
    hptframe->SetXTitle("E_{T} (GeV/c)"); 
    hptframe->SetYTitle("#frac{1}{N_{evt}}#frac{d^{2}N}{dE_{T}d#eta }");
  }
  hptframe->Draw();
  //  TH1D * pt[12];
  //  TGraphErrors * gpt[12];
  //  TLegend * leg2 = new TLegend(0.60,0.50,0.85,0.88,"Centrality    (N_{part})");
  for(int icent = 0; icent<MaxCentBin; icent++) {
//     if(type==1) {
//       pt[icent] = (TH1D *) dNdPt->ProjectionX(Form("hpt_%d_%d",mincent[icent],maxcent[icent]),icent+1,icent+1);
//     } else { 
//       pt[icent] = (TH1D *) dNdEt->ProjectionX(Form("het_%d_%d",mincent[icent],maxcent[icent]),icent+1,icent+1);
//     }
//     int nptbins = pt[icent]->GetNbinsX();
//     int ncnt = 0;
//     for(int i = 0; i<nptbins; i++ ) {
//       if(pt[icent]->GetBinContent(i+1)>0) {
// 	Double_t scale = 1./hcent->Integral(hcent->FindBin(mincent[icent]),hcent->FindBin(maxcent[icent]-0.1));
// 	yy[ncnt]= scale*pt[icent]->GetBinContent(i+1);
// 	yyerr[ncnt]=scale*pt[icent]->GetBinError(i+1);
// 	xxerr[ncnt]=0;
// 	if(type==1) {
// 	  xx[ncnt]=hpt->GetBinContent(i+1,icent+1);
// 	} else if (type==2 || type==3) {
// 	  xx[ncnt]=het->GetBinContent(i+1,icent+1);
//       }
// 	++ncnt;
//       } 
      
//     }
//     gpt[icent] = new TGraphErrors(ncnt,xx,yy,xxerr,yyerr);
//     gpt[icent]->SetMarkerStyle(markers[icent]);
//     gpt[icent]->SetMarkerColor(colors[icent]);
//     gpt[icent]->SetLineColor(colors[icent]);
   gpt[icent]->Draw("p");
   //leg2->AddEntry(g[icent],  Form("%d-%d     (%5.1f)",mincent[icent],maxcent[icent],NpartBin->GetBinContent(icent+1)),"lp");
  }
  leg2->Draw();   
  TPaveText * desc2;
  if(type == 1 ) {
    desc2 = new TPaveText(0.2,0.01,2.8,0.2);
  } else if (type == 2||type==3) {
    desc2 = new TPaveText(0.2,0.01,2.8,0.2);
  }

  desc2->AddText(label.Data());

  if(type==1) {
    desc2->AddText(Form("%5.1f #leq #eta_{track} < %5.1f",mineta1,maxeta1));
  } else if (type==2||type==3) {
    desc2->AddText(Form("%5.1f #leq #eta_{calo} < %5.1f",mineta1,maxeta1));
  }
  if(mineta2<maxeta2) {
    desc2->AddText(Form("RP: %s",EPNames[EP2].data()));
    if(type==1) {
      desc2->AddText(Form("%5.1f #leq #eta_{track} < %5.1f",mineta2,maxeta2));
    } else if (type==2||type==3) {
      desc2->AddText(Form("%5.1f #leq #eta_{calo} < %5.1f",mineta2,maxeta2));
    }
  }
  desc2->Draw();

  if(mineta2==maxeta2) {
    c3->Print(Form("~/public_html/PtSpectra/ptDist_%s_%s_%s_%s.png",typeName.Data(),EPNames[EP1].data(),EtaBins.Data(),tag.Data()),"png");
  } else {
    c3->Print(Form("~/public_html/PtSpectra/ptDist_%s_%s_%s_%s_%s.png",typeName.Data(),EPNames[EP1].data(),EPNames[EP2].data(),EtaBins.Data(),tag.Data()),"png");
  }

  TString intV2FileName = "results/intv2_";
  if(effFileName.Contains("bass20k")) intV2FileName+="bass20k_";
  if(effFileName.Contains("ampt30k")) intV2FileName+="ampt30k_";
  if(effFileName.Contains("NoEffCorr")) intV2FileName+="NoEffCorr_";
  intV2FileName+=Form("%03d_%03d_",(int)(10.*mineta1),(int)(10.*maxeta1));
  if(maxeta2>mineta2) intV2FileName+=Form("%03d_%03d_",(int)(10.*mineta2),(int)(10.*maxeta2));
  intV2FileName+=tag.Data();
  intV2FileName+=".txt";
  std::ofstream fileInt;
  fileInt.open(intV2FileName.Data());
  double err=0;
  double dndetaInt = 0;
  double dndetaerrInt=0;
  fout->cd();
  if(MINETA2==MAXETA2) {
    fout->mkdir(Form("%04.1f_%04.1f",MINETA1,MAXETA1));
  } else {
    fout->mkdir(Form("%04.1f_%04.1f_%04.1f_%04.1f",MINETA1,MAXETA1,MINETA2,MAXETA2));
  }
  for(int icent = 0; icent<MaxCentBin; icent++) {
    MINCENT = mincent[icent];
    MAXCENT = maxcent[icent];
    double v2int = IntegralV2(pt[icent], g[icent],gpt[icent],minpt,maxpt,err,dndetaInt,dndetaerrInt);
    cout<<"Integral v2: "<<mincent[icent]<<"-"<<maxcent[icent]<<" = "<<v2int<<" +/- "<<err<<"  dNdEta: "<<dndetaInt<<endl;
    fileInt<<setprecision(3)<<aveta1<<"\t"<<setprecision(3)<<v2int<<"\t"<<setprecision(3)<<err<<"\t"<<setprecision(3)<<mincent[icent]<<"\t"<<setprecision(3)<<maxcent[icent]<<"\t"<<setprecision(4)<<dndetaInt<<"\t"<<setprecision(3)<<dndetaerrInt<<"\t"<<endl;
  }
  fileInt.close();
}

double IntegralV2(TH1D * pt,TGraphErrors * g,TGraphErrors * gpt, double ptmin, double ptmax, double &err, double &dndeta, double &dndetaerr){
  double ret = 0;
  double norm = 0;
  
  double err2 = 0;
  double err2pt = 0;
  double * gx;
  double * gy;
  double * gey;
  double * gptx;
  double * gpty;
  double * gpte;
  TString specname;
  TString v2name;
  TString gspecname;
  TString gv2name;
  fout->cd();
  if(MINETA2==MAXETA2) {
    fout->cd( Form("%04.1f_%04.1f",MINETA1,MAXETA1));
    specname=Form("hdNdPt_%d_%d",(Int_t)MINCENT,(Int_t)MAXCENT);
    v2name=  Form("hv2_%d_%d",     (Int_t)MINCENT,(Int_t)MAXCENT);
    gspecname=Form("dNdPt_%d_%d",(Int_t)MINCENT,(Int_t)MAXCENT);
    gv2name=  Form("v2_%d_%d",      (Int_t)MINCENT,(Int_t)MAXCENT);
  } else {
    fout->cd( Form("%04.1f_%04.1f_%04.1f_%04.1f",MINETA1,MAXETA1,MINETA2,MAXETA2));
    specname=Form("hdNddPt_%d_%d",(Int_t)MINCENT,(Int_t)MAXCENT);
    v2name=Form("hv2_%d_%d",       (Int_t)MINCENT,(Int_t)MAXCENT);
    gspecname=Form("dNdPt_%d_%d",(Int_t)MINCENT,(Int_t)MAXCENT);
    gv2name=Form("v2_%d_%d",       (Int_t)MINCENT,(Int_t)MAXCENT);
  }
  gpt->SetName(gspecname.Data());
  g->SetName(gv2name.Data());
  gx = g->GetX();
  gy = g->GetY();
  gey = g->GetEY();
  TH1D * hspec = new TH1D(specname.Data(),specname.Data(),100,0,12);
  hspec->SetYTitle("#frac{1}{N_{evt}}#frac{d^{2}N}{dp_{T}d#eta }");
  hspec->SetXTitle("p_{T} (GeV/c)");
  TH1D * hv2 = new TH1D(v2name.Data(),v2name.Data(),100,0,12);
  hv2->SetXTitle("p_{T} (GeV/c)");
  hv2->SetYTitle("v_{2}(EP)");
  int gptn = g->GetN();
  if(ptmin>0.2) {
    for(int i = 0; i<g->GetN(); i++) {
      gx[i]=gx[i+1];
      gy[i]=gy[i+1];
      gey[i]=gey[i+1];
    }
    gptn-=1;
  }
  TGraphErrors * newg = new TGraphErrors(gptn,gx,gy,0,gey);
  newg->SetName(gv2name.Data());
  newg->SetTitle(gv2name.Data());
  newg->Write();
  newg->SetHistogram(hv2);
  newg->GetHistogram()->SetMinimum(0.0);
  newg->GetHistogram()->SetMaximum(0.4);
  newg->GetHistogram()->SetXTitle("p_{T} (GeV/c)");
  newg->GetHistogram()->SetYTitle("v_{2}(EP)");
  gptx = gpt->GetX();
  gpty = gpt->GetY();
  gpte = gpt->GetEY();
  gptn = gpt->GetN();
  if(ptmin>0.2) {
    for(int i = 0; i<gpt->GetN(); i++) {
      gptx[i]=gptx[i+1];
      gpty[i]=gpty[i+1];
      gpte[i]=gpte[i+1];
    }
    gptn -= 1;
  }
  TGraphErrors * newgpt = new TGraphErrors(gptn,gptx,gpty,0,gpte);
  newgpt->SetName(gspecname.Data());
  newgpt->SetTitle(gspecname.Data());
  newgpt->Write();
  hspec->SetMinimum(0.01);
  hspec->SetMaximum(10000.);
  newgpt->SetHistogram(hspec);
  err2 = 0;
  for(Int_t i = 0; i<newg->GetN(); i++) {
    int ptbin = pt->FindBin(gx[i]);
    double dpt = pt->GetBinWidth(ptbin);
    if(pt->GetBinLowEdge(ptbin)>=ptmin && pt->GetBinLowEdge(ptbin)<ptmax) {
      ret +=      gy[i]*gpty[i]*dpt;
      err2 +=pow(gey[i]*gpty[i]*dpt,2);
      err2pt+=pow(gpte[i]*dpt,2);
      norm+=gpty[i]*dpt;
    }
  }

  tf->cd();
  err=sqrt(err2)/norm;
  ret/=norm;
  dndeta = norm;
  dndetaerr=sqrt(err2pt);
  return ret;
}
void AddToV2(Int_t type, Int_t rp, int minc, int maxc, double mineta, double maxeta, TH1D ** rescor) {
  TString name = EPNames[rp];
  Double_t mincent = minc+0.01;
  Double_t maxcent = maxc-0.01;
  Int_t icentmin = hcentbins->FindBin(mincent);
  Int_t icentmax = hcentbins->FindBin(maxcent-1.);
  Int_t ietamin =  hetabins->FindBin(mineta);
  Int_t ietamax = hetabins->FindBin(maxeta-0.1);
  TString cosname; 
  TString cntname; 
  TString prefix="";
  if(GenCalc) prefix="gen";
  Int_t ioff = 0;
  if(!cos_) {
    if(type ==1 ) {
      cosname = Form("v2analyzer/v2/v2Reco/%s/%scos_%s_%d",name.Data(),prefix.Data(),name.Data(),ietamin-1);
      cntname = Form("v2analyzer/v2/v2Reco/%s/%scnt_%s_%d",name.Data(),prefix.Data(),name.Data(),ietamin-1);
    } else if (type == 2||type==3) {
      cosname = Form("v2analyzer/v2/v2Reco/calo_%s/%scos_calo_%s_%d",name.Data(),prefix.Data(),name.Data(),ietamin-1);
      cntname = Form("v2analyzer/v2/v2Reco/calo_%s/%scnt_calo_%s_%d",name.Data(),prefix.Data(),name.Data(),ietamin-1);
    }
    cos_ = (TH2D *)((TH2D *) tf->Get(cosname.Data()))->Clone(Form("cos_%s_%d_%d",name.Data(),icentmin,ietamin));
    cnt_ = (TH2D *)((TH2D *) tf->Get(cntname.Data()))->Clone(Form("cnt_%s_%d_%d",name.Data(),icentmin,ietamin));
    ioff = 1;
  }
  for(Int_t ib = ietamin+ioff; ib<=ietamax; ib++) {
    if(type == 1 ) {
      cos_->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/%s/%scos_%s_%d",name.Data(),prefix.Data(),name.Data(),ib-1)));
      cnt_->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/%s/%scnt_%s_%d",name.Data(),prefix.Data(),name.Data(),ib-1)));
    } else if(type==2||type==3) {
      cos_->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/calo_%s/%scos_calo_%s_%d",name.Data(),prefix.Data(),name.Data(),ib-1)));
      cnt_->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/calo_%s/%scnt_calo_%s_%d",name.Data(),prefix.Data(),name.Data(),ib-1)));
    }
  }
  for(Int_t i = 1; i<= cos_->GetNbinsX(); i++) {
    Double_t resc = rescor[rp]->GetBinContent(i);
    if(GenCalc) resc = 1.;
    for(Int_t j = 1; j<=cos_->GetNbinsY();j++) {
      if(resc>0) {
	cos_->SetBinContent(i,j,cos_->GetBinContent(i,j)/resc);
	cos_->SetBinError(i,j,cos_->GetBinError(i,j)/resc);
      } else {
	cos_->SetBinContent(i,j,0.);
      }
    }
  } 
}
TGraphErrors * GenV2(double mincent, double maxcent, TH2D * hpt, Int_t marker, Int_t color) {
  Int_t icentmin = hcentbins->FindBin(mincent);
  Int_t icentmax = hcentbins->FindBin(maxcent-1.);
  cos_->Divide(cnt_);
  TH1D * v2pt = (TH1D *) cos_->ProjectionY(Form("v2pt_%d_%d",icentmin,icentmax),icentmin,icentmax);
  v2pt->Reset();
  for(int i = 1; i<=v2pt->GetNbinsX(); i++) { 
    double err = 0;
    for(int j = icentmin; j<=icentmax; j++) {
      v2pt->SetBinContent(i,v2pt->GetBinContent(i)+cos_->GetBinContent(j,i));
      err+=cos_->GetBinError(j,i)*cos_->GetBinError(j,i);
    }
    if(err>0) err = sqrt(err);
    v2pt->SetBinError(i,err);
  }
  Double_t xx[nPtBins];
  Double_t xxerr[nPtBins];
  Double_t yy[nPtBins];
  Double_t yyerr[nPtBins];
  Int_t npt=0;
  TH1D * hcent;
  hptcent = (TH1D *) hpt->ProjectionX(Form("tmp%s_%d_%d",hpt->GetName(),icentmin,icentmax),icentmin,icentmax);
  for(int i = 0; i< v2pt->GetNbinsX(); i++ ) {
    if(hptcent->GetBinContent(i+1)>0) {
      yy[npt]    = v2pt->GetBinContent(i+1);
      yyerr[npt] = v2pt->GetBinError(i+1);
      xx[npt]    = hpt->GetBinContent(i+1,icentmin);
      xxerr[npt] = 0;
      //if(xx[npt]>minpt) ++npt;
      ++npt;
    }
  }
  TGraphErrors * GenV2 = new TGraphErrors(npt,xx,yy,xxerr,yyerr);
  GenV2->SetMarkerStyle(marker);
  GenV2->SetMarkerColor(color);
  GenV2->SetLineColor(color);
  return GenV2;
}
