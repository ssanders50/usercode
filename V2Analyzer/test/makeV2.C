#include "/net/hisrv0001/home/sanders/CMSSW_3_9_1/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include <iomanip>
#include <fstream>
using namespace std;
TFile * tf;
double centbins[]={0,5,10,15,20,30,40,50,60,70,80,90,100};
Int_t nCentBins = 12;
double ptbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0};
Int_t nPtBins = 15;
double etabins[]={-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5};
Int_t nEtaBins = 20;
TH1D * hcentbins = new TH1D("centbins","centbins",nCentBins,centbins);
TH1D * hptbins = new TH1D("ptbins","ptbins",nPtBins,ptbins);
TH1D * hetabins = new TH1D("hetabins","hetabins",nEtaBins,etabins);

TGraphErrors * GenV2(Int_t type, Int_t rp, int minc, int maxc,  Double_t mineta, Double_t maxeta, TH1D ** rescor, TH2D * hpt, Int_t marker, Int_t color);
Bool_t GenCalc = kFALSE;
void makeV2(){
  //--------------------------------------
  // Set up analysis
  //
  //  TString tag = "Hydjet_Bass";
  TString tag = "Hydjet_Bass";
  if(GenCalc) tag = tag+"_GENERATOR";
  Int_t EP = etHFp;
  //Int_t EP = EvtPTracksPosEtaGap;
  //Int_t EP = EvtPlaneFromTracksMidEta;
  //Int_t EP = etCaloHFP;
  Int_t type = 1;  // =1 for tracks, =2 for etcalo, = 3 for hcal(special replay)
  double mineta = -2.;
  double maxeta = -1. ;
  //
  // End of setup
  //-------------------------------------
  Int_t ietamin =  hetabins->FindBin(mineta);
  Int_t ietamax = hetabins->FindBin(maxeta-0.1);
  TString EtaBins = Form("etabins%d_%d",ietamin,ietamax);
  mineta = hetabins->GetBinLowEdge(ietamin);
  maxeta = hetabins->GetBinLowEdge(ietamax)+hetabins->GetBinWidth(ietamax);
  //                0-5   5-10 10-15 15-20   20-30   30-40   40-50   50-60  60-70    70-80 80-90 90-100
  int mincent[12] ={0,     5,    10,  15,     20,     30,     40,     50,    60,      70,   80,    90       };
  int maxcent[12] ={5,    10,    15,  20,     30,     40,     50,     60,    70,      80,   90,   100       };
  int markers[12]= {22,    21,   23,  24,     25,     20,     29,     28,    26,      27,   3,     5        };
  int colors[12] = {kBlack,kBlue,kRed,kCyan+2,kViolet,kSpring,kOrange,kRed+4,kAzure+9,kTeal,kRed-4,kYellow+2};

  tf = new TFile("/net/hisrv0001/home/sanders/CMSSW_3_9_1/src/V2Analyzer/V2Analyzer/data/rpflat_combined.root");

  // Step 0.  Generate Npart is spectra
  
  TH2D * hpt = (TH2D *) tf->Get(Form("v2analyzer/pt_%d",ietamin-1))->Clone("hpt");
  TH2D * dNdPt = (TH2D *) tf->Get(Form("v2analyzer/ptCnt_%d",ietamin-1))->Clone("dNdPt");
  TH2D * het = (TH2D *) tf->Get(Form("v2analyzer/et_%d",ietamin-1))->Clone("hpt");
  TH2D * dNdEt = (TH2D *) tf->Get(Form("v2analyzer/etCnt_%d",ietamin-1))->Clone("dNdPt");
  for(int ieta = ietamin; ieta<ietamax; ieta++) {
    hpt->Add((TH2D *) tf->Get(Form("v2analyzer/pt_%d",ieta)));
    dNdPt->Add((TH2D *) tf->Get(Form("v2analyzer/ptCnt_%d",ieta)));
    het->Add((TH2D *) tf->Get(Form("v2analyzer/et_%d",ieta)));
    dNdEt->Add((TH2D *) tf->Get(Form("v2analyzer/etCnt_%d",ieta)));
  }
  hpt->Divide(dNdPt);
  for(int i = 1; i<= dNdPt->GetNbinsX(); i++ ) 
    for(int j = 1; j<=dNdPt->GetNbinsY(); j++) {
      dNdPt->SetBinContent(i,j, dNdPt->GetBinContent(i,j)/
			   dNdPt->GetXaxis()->GetBinWidth(i)/2./TMath::Pi());
    }

  het->Divide(dNdEt);
  for(int i = 1; i<= dNdEt->GetNbinsX(); i++ ) 
    for(int j = 1; j<=dNdEt->GetNbinsY(); j++) {
      dNdEt->SetBinContent(i,j, dNdEt->GetBinContent(i,j)/
			   dNdEt->GetXaxis()->GetBinWidth(i)/2./TMath::Pi());
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
  Double_t xx[25];
  Double_t yy[25];
  Double_t xxerr[25];
  Double_t yyerr[25];
  TCanvas * c2 = new TCanvas("c2","c2",800,600);
  gPad->SetGrid(1,1);

  TH1D * hfig;
  if(type == 1 ) {
    hfig = new TH1D("hfig",Form("%s",tag.Data()),100,0,6);
  } else if (type == 2 || type == 3) {
    hfig = new TH1D("hfig",Form("%s",tag.Data()),100,0,12);
  }
  hfig->SetMinimum(-0.05);
  hfig->SetMaximum(0.30);
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
    desc = new TPaveText(0.1,0.25,2.6,0.29);
  } else if (type == 2 || type==3) {
    desc = new TPaveText(0.1,0.25,3.8,0.29);
  }
  TString label = Form("RP: %s",EPNames[EP]);
  if(GenCalc) label = "GENERATOR";
  desc->AddText(label.Data());
  if(type==1) {
    desc->AddText(Form("%5.1f #leq #eta_{track} < %5.1f",mineta,maxeta));
  } else if (type==2) {
      desc->AddText(Form("%5.1f #leq #eta_{calo} < %5.1f",mineta,maxeta));
  } else if (type==3) {
      desc->AddText(Form("%5.1f #leq #eta_{hcal} < %5.1f",mineta,maxeta));
  }
  desc->Draw();
  TGraphErrors * g[12];
  TLegend * leg = new TLegend(0.65,0.6,0.85,0.88,"Centrality    (N_{part})");
  std::ofstream file;
  file.open(Form("results/v2Results_%s_%s_%s_%s.txt",typeName.Data(),EPNames[EP],EtaBins.Data(),tag.Data()));
  file<<tag.Data()<<endl<<endl;
  for(Int_t i = 0; i<8; i++) {
    if(type == 1 ) {
      g[i]   = GenV2(type,EP,mincent[i],maxcent[i], mineta, maxeta, rescor,hpt, markers[i],  colors[i]);
    } else if (type == 2||type==3) {
      g[i]   = GenV2(type,EP,mincent[i],maxcent[i], mineta, maxeta, rescor,het, markers[i],  colors[i]);
    }
    g[i]->Draw("p");
    leg->AddEntry(g[i],  Form( "%d-%d     (%5.1f)",mincent[i],maxcent[i],NpartBin->GetBinContent(i+1)),"lp");
    file<<Form("%d-%d     (Npart = %5.1f)",mincent[i],maxcent[i],NpartBin->GetBinContent(1))<<endl;
    Double_t * xxx = g[i]->GetX();
    Double_t * yyy = g[i]->GetY();  
    Double_t * yyyerr = g[i]->GetEY();
        for(int j = 0; j< g[i]->GetN();  j++ ) {
      if(yyy[j]<0.0001) continue;
      file<<setprecision(3)<<xxx[j]<<"\t"<<setprecision(3)<<yyy[j]<<"\t"<<yyyerr[j]<<endl;
     }
    file<<endl;
  } 
  leg->Draw();
  c2->Print(Form("~/public_html/V2Spectra/v2_%s_%s_%s_%s.png",typeName.Data(),EPNames[EP].data(),EtaBins.Data(),tag.Data()),"png");
  
  TCanvas * c3; 
  if(type==1 ) {
    c3= new TCanvas("ptDist","ptDist",800,600);
  } else if (type==2||type==3) {
    c3= new TCanvas("etDist","etDist",800,600);
  }
  TH1D * hptframe;
  if(type == 1 ) {
    hptframe = new TH1D("hptframe",tag.Data(),100,0,6);
  } else if (type == 2 || type==3) {
    hptframe = new TH1D("hptframe",tag.Data(),100,0,12);
  }
  hptframe->SetMaximum(10000000);
  hptframe->SetMinimum(0.01);
  hptframe->SetStats(kFALSE);
  gPad->SetLogy();
  if(type==1) {
    hptframe->SetXTitle("p_{T} (GeV/c)"); 
    hptframe->SetYTitle("#frac{1}{2#pi p_{t}}#frac{d^{2}N}{dp_{T}d#eta }");
  } else if(type==2||type==3) {
    hptframe->SetXTitle("E_{T} (GeV/c)"); 
    hptframe->SetYTitle("#frac{1}{2#pi E_{t}}#frac{d^{2}N}{dE_{T}d#eta }");
  }
  hptframe->Draw();
  TH1D * pt[12];
  TGraphErrors * gpt[12];
  TLegend * leg2 = new TLegend(0.60,0.50,0.85,0.88,"Centrality    (N_{part})");
  for(int icent = 0; icent<8; icent++) {
    pt[icent] = (TH1D *) dNdPt->ProjectionX(Form("hpt_%d_%d",mincent[icent],maxcent[icent]),1,1);
    int nptbins = pt[icent]->GetNbinsX();
    for(int i = 0; i<nptbins; i++ ) {
      xx[i]=hpt->GetBinContent(i+1,hcentbins->FindBin(mincent[icent]+0.1),hcentbins->FindBin(maxcent[icent]-0.1));
      if(xx[i]>0) {
	Double_t scale = 1./(2.*3.1415*xx[i]);
	yy[i]= scale*pt[icent]->GetBinContent(i+1)/pow(2,icent);
	yyerr[i]=scale*pt[icent]->GetBinError(i+1);
      } else {
	yy[i] = 0;
	yyerr[i] = 0;
      }
      
    }
    gpt[icent] = new TGraphErrors(nptbins,xx,yy,xxerr,yyerr);
    gpt[icent]->SetMarkerStyle(markers[icent]);
    gpt[icent]->SetMarkerColor(colors[icent]);
    gpt[icent]->SetLineColor(colors[icent]);
    gpt[icent]->Draw("p");
    leg2->AddEntry(g[icent],  Form("%d-%d     (%5.1f)",mincent[icent],maxcent[icent],NpartBin->GetBinContent(1)),"lp");
  }
  leg2->Draw();   
  TPaveText * desc2;
  if(type == 1 ) {
    desc2 = new TPaveText(0.2,0.1,1.8,5.0);
  } else if (type == 2||type==3) {
    desc2 = new TPaveText(0.2,0.1,2.8,5.0);
  }
  desc2->AddText(Form("P: %s",EPNames[EP]));
  if(type==1) {
    desc2->AddText(Form("%5.1f #leq #eta_{track} < %5.1f",mineta,maxeta));
  } else if (type==2||type==3) {
    desc2->AddText(Form("%5.1f #leq #eta_{calo} < %5.1f",mineta,maxeta));
  }
  desc2->Draw();


  c3->Print(Form("~/public_html/PtSpectra/ptDist_%s_%s_%s_%s.png",typeName.Data(),EPNames[EP].data(),EtaBins.Data(),tag.Data()),"png");
  file.close();
}

TGraphErrors * GenV2(Int_t type, Int_t rp, int minc, int maxc,  Double_t mineta, Double_t maxeta, TH1D ** rescor, TH2D * hpt, Int_t marker, Int_t color) {
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
  if(type ==1 ) {
    cosname = Form("v2analyzer/v2/v2Reco/%s/%scos_%s_%d",name.Data(),prefix.Data(),name.Data(),ietamin-1);
    cntname = Form("v2analyzer/v2/v2Reco/%s/%scnt_%s_%d",name.Data(),prefix.Data(),name.Data(),ietamin-1);
  } else if (type == 2||type==3) {
    cosname = Form("v2analyzer/v2/v2Reco/calo_%s/%scos_calo_%s_%d",name.Data(),prefix.Data(),name.Data(),ietamin-1);
    cntname = Form("v2analyzer/v2/v2Reco/calo_%s/%scnt_calo_%s_%d",name.Data(),prefix.Data(),name.Data(),ietamin-1);
  }
  TH2D * cos_ = (TH2D *)((TH2D *) tf->Get(cosname.Data()))->Clone(Form("cos_%s_%d_%d",name.Data(),icentmin,ietamin));
  TH2D * cnt_ = (TH2D *)((TH2D *) tf->Get(cntname.Data()))->Clone(Form("cnt_%s_%d_%d",name.Data(),icentmin,ietamin));
  for(Int_t ib = ietamin+1; ib<=ietamax; ib++) {
    if(type == 1 ) {
      cos_->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/%s/%scos_%s_%d",name.Data(),prefix.Data(),name.Data(),ib-1)));
      cnt_->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/%s/%scnt_%s_%d",name.Data(),prefix.Data(),name.Data(),ib-1)));
    } else if(type==2++type==3) {
      cos_->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/calo_%s/%scos_calo_%s_%d",name.Data(),prefix.Data(),name.Data(),ib-1)));
      cnt_->Add((TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/calo_%s/%scnt_calo_%s_%d",name.Data(),prefix.Data(),name.Data(),ib-1)));
    }
  }
  cos_->Divide(cnt_);
  for(Int_t i = 1; i<= cos_->GetNbinsX(); i++) {
    Double_t resc = rescor[rp]->GetBinContent(i);
    if(GenCalc) resc = 1.;
    for(Int_t j = 1; j<=cos_->GetNbinsY();j++) {
      if(resc>0) {
	cos_->SetBinContent(i,j,cos_->GetBinContent(i,j)/resc);
      } else {
	cos_->SetBinContent(i,j,0.);
      }
    }
  }
  TH1D * v2pt = (TH1D *) cos_->ProjectionY(Form("v2pt_%s_%d_%d",name.Data(),icentmin,icentmax),icentmin,icentmax);
  Double_t xx[40];
  Double_t xxerr[40];
  Double_t yy[40];
  Double_t yyerr[40];
  Int_t npt=0;
  for(int i = 0; i< v2pt->GetNbinsX(); i++ ) {
   
    if(v2pt->GetBinContent(i+1) > 0&&hpt->GetBinContent(i+1,icentmin)>0) {
      yy[npt]    = v2pt->GetBinContent(i+1);
      yyerr[npt] = v2pt->GetBinError(i+1);
      xx[npt]    = hpt->GetBinContent(i+1,icentmin);
      xxerr[npt] = 0;
      ++npt;
    }
  }
  TGraphErrors * GenV2 = new TGraphErrors(npt,xx,yy,xxerr,yyerr);
  GenV2->SetMarkerStyle(marker);
  GenV2->SetMarkerColor(color);
  GenV2->SetLineColor(color);
  return GenV2;
}
