#include "/net/hisrv0001/home/sanders/CMSSW_3_9_2_patch5/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include <iomanip>
#include <fstream>
using namespace std;
TFile * tf;

static const Int_t nCentBins = 12;
static const Int_t nPtBins = 15;
static const Int_t nEtBins = 15;
static const Int_t nEtaBins = 22;
static const double minpt = 0.3;
static const double centbins[]={0,5,10,15,20,30,40,50,60,70,80,90,100};
static const double ptbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0};
static const double etbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0};
static const double etabins[]={-5,-4.5,-4,-3.5,-3,-2.4,-2,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,3,3.5,4,4.5,5};
// STANDARD ETA BINS:     |eta|<0.8
//                  0.8 < |eta| < 1.6
//                  1.6 < |eta| < 2.0
//                  2.0 < |eta| < 2.4
TH1D * hcentbins = new TH1D("centbins","centbins",nCentBins,centbins);
TH1D * hptbins = new TH1D("ptbins","ptbins",nPtBins,ptbins);
TH1D * hetabins = new TH1D("hetabins","hetabins",nEtaBins,etabins);

void AddToV2(Int_t type, Int_t rp, int minc, int maxc,  Double_t mineta, Double_t maxeta, TH1D ** rescor);
double IntegralV2(Int_t icent,TH1D * pt,TGraphErrors * g,TGraphErrors * gpt, double &err);
TGraphErrors * GenV2(double centmin, double centmax, TH2D * hpt, Int_t marker, Int_t color);
Bool_t GenCalc = kFALSE;
TH2D * cos1_;
TH2D * cnt1_;
TH2D * cos_;
TH2D * cnt_;
TFile * feff;
void makeV2(){
  //--------------------------------------
  // Set up analysis
  //
  if(GenCalc) tag = tag+"_GENERATOR";
  Int_t type = 1;  // =1 for tracks, =2 for etcalo, = 3 for hcal(special replay)
  Int_t EP1 = EPTracksPosEtaBigGap;
  double mineta1 = -0.8;
  double maxeta1 = 0.0;
  double aveta1 = (fabs(mineta1)+fabs(maxeta1))/2.;
  Int_t EP2 = EPTracksNegEtaBigGap;
  double mineta2 = 0.0;
  double maxeta2 = 0.8;
  TString tag = "_hiGoodMergedTracks_03_25PtCut_PtWeight";
  Int_t MaxCentBin = 9;
  //
  // End of setup
  //-------------------------------------
  cos_ = 0;
  cnt_ = 0;
  cos1_ = 0;
  cnt1_ = 0;
  Int_t ietamin1 =  hetabins->FindBin(mineta1);
  Int_t ietamax1 = hetabins->FindBin(maxeta1-0.1);
  Int_t ietamin2 = hetabins->FindBin(mineta2);
  Int_t ietamax2 = hetabins->FindBin(maxeta2-0.1);
  TString EtaBins = Form("etabins_%d_%d",ietamin1,ietamax1);
  mineta1 = hetabins->GetBinLowEdge(ietamin1);
  maxeta1 = hetabins->GetBinLowEdge(ietamax1)+hetabins->GetBinWidth(ietamax1);
  if(mineta2<maxeta2) {
    EtaBins += Form("_%d_%d",ietamin2,ietamax2);
    mineta2 = hetabins->GetBinLowEdge(ietamin2);
    maxeta2 = hetabins->GetBinLowEdge(ietamax2)+hetabins->GetBinWidth(ietamax2);
  }

  cout<<"Eta ranges: "<<mineta1<<"-"<<maxeta1<<" :  "<<mineta2<<"-"<<maxeta2<<endl;
  //                0-5   5-10 10-15 15-20   20-30   30-40   40-50   50-60  60-70    70-80 80-90 90-100
  int mincent[12] ={0,     5,    10,  15,     20,     30,     40,     50,    60,      70,   80,    90       };
  int maxcent[12] ={5,    10,    15,  20,     30,     40,     50,     60,    70,      80,   90,   100       };
  int markers[12]= {22,    21,   23,  24,     25,     20,     29,     28,    26,      27,   3,     5        };
  int colors[12] = {kBlack,kBlue,kRed,kCyan+2,kViolet,kSpring,kOrange,kRed+4,kAzure+9,kTeal,kRed-4,kYellow+2};

  tf = new TFile("$CMSSW_BASE/src/V2Analyzer/V2Analyzer/data/rpflat_combined.root");

  // Step 0.  Generate Npart is spectra
  
  TH2D * hpt = (TH2D *) tf->Get(Form("v2analyzer/Spectra/pt_%d",ietamin1-1))->Clone("hpt");
  TH2D * dNdPt = (TH2D *) tf->Get(Form("v2analyzer/Spectra/ptCnt_%d",ietamin1-1))->Clone("dNdPt");
  TH2D * het = (TH2D *) tf->Get(Form("v2analyzer/Spectra/et_%d",ietamin1-1))->Clone("het");
  TH2D * dNdEt = (TH2D *) tf->Get(Form("v2analyzer/Spectra/etCnt_%d",ietamin1-1))->Clone("dNdEt");
  TH1D * hcent = (TH1D *) tf->Get("v2analyzer/cent");
  double deta  = etabins[ietamin1]-etabins[ietamin1-1];
  double deta2 = 0;
  cout<<"Add: "<<ietamin1-1<<endl;
  TH2D * heff1 = (TH2D *) dNdPt->Clone("heff");
  TH2D * heff2 = (TH2D *) dNdPt->Clone("heff");
  heff1->Reset();
  heff2->Reset();
  feff = new TFile("trkCorrFlow_HYDJET_20k.root");

  //Note..the following assumes request for the "standard" eta ranges, which is all Eric supplies in his efficiency maps.
  for(int icent = 0; icent< MaxCentBin; icent++) {
    TH2D * heffRef = (TH2D *) feff->Get(Form("rEff_cbin%d",icent));
    int ieta1 = heffRef->GetXaxis()->FindBin(mineta1+0.1);
    int ieta2 = 0;
    if(mineta2!=maxeta2) ieta2=heffRef->GetXaxis()->FindBin(mineta2+0.1);
    for(int ipt = 0; ipt<14; ipt++) {
      heff1->SetBinContent(ipt+2,icent+1,heffRef->GetBinContent(ieta1,ipt+1));
      if(ieta2>0) heff2->SetBinContent(ipt+2,icent+1,heffRef->GetBinContent(ieta2,ipt+1));
    }
  }
  if(ieta2>0) {
    heff1->Add(heff2);
    heff1->Scale(0.5);
  }
  for(int ieta = ietamin1; ieta<ietamax1; ieta++) {
    deta+=etabins[ieta+1]-etabins[ieta];
    hpt->Add((TH2D *) tf->Get(Form("v2analyzer/Spectra/pt_%d",ieta)));
    dNdPt->Add((TH2D *) tf->Get(Form("v2analyzer/Spectra/ptCnt_%d",ieta)));
    het->Add((TH2D *) tf->Get(Form("v2analyzer/Spectra/et_%d",ieta)));
    dNdEt->Add((TH2D *) tf->Get(Form("v2analyzer/Spectra/etCnt_%d",ieta)));
    cout<<"Add: "<<ieta<<endl;
  }
  if(mineta2<maxeta2) {
    for(int ieta = ietamin2-1; ieta<ietamax2; ieta++) {
      deta2+=etabins[ieta+1]-etabins[ieta];
      hpt->Add((TH2D *) tf->Get(Form("v2analyzer/Spectra/pt_%d",ieta)));
      dNdPt->Add((TH2D *) tf->Get(Form("v2analyzer/Spectra/ptCnt_%d",ieta)));
      het->Add((TH2D *) tf->Get(Form("v2analyzer/Spectra/et_%d",ieta)));
      dNdEt->Add((TH2D *) tf->Get(Form("v2analyzer/Spectra/etCnt_%d",ieta)));
      cout<<"Add: "<<ieta<<endl;
    }
  }
  if(deta!=deta2) cout<<"WARNING:  ETA Assymmetry! with deta/deta2= "<<deta<<"/"<<deta2<<endl;
  hpt->Divide(dNdPt);
  for(int i = 1; i<= dNdPt->GetNbinsX(); i++ ) 
    for(int j = 1; j<=dNdPt->GetNbinsY(); j++) {
      dNdPt->SetBinContent(i,j, dNdPt->GetBinContent(i,j)/
			   (dNdPt->GetXaxis()->GetBinWidth(i)*deta));
    }
  dNdPt->Divide(heff1);
  het->Divide(dNdEt);
  for(int i = 1; i<= dNdEt->GetNbinsX(); i++ ) 
    for(int j = 1; j<=dNdEt->GetNbinsY(); j++) {
      dNdEt->SetBinContent(i,j, dNdEt->GetBinContent(i,j)/
			   (dNdEt->GetXaxis()->GetBinWidth(i)*deta));
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
    hfig = new TH1D("hfig",Form("%s",tag.Data()),100,0,9);
  } else if (type == 2 || type == 3) {
    hfig = new TH1D("hfig",Form("%s",tag.Data()),100,0,12);
  }
  hfig->SetMinimum(0.0);
  if(type==2) hfig->SetMinimum(-0.1);
  hfig->SetMaximum(0.3);
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
  if(GenCalc) label = "GENERATOR";
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
  TGraphErrors * g[12];
  TLegend * leg = new TLegend(0.75,0.67,0.92,0.92,"Centrality    (N_{part})");
  std::ofstream file;
  double emax = fabs(mineta1);
  if(fabs(maxeta1)>emax) emax = fabs(maxeta1);
  for(Int_t i = 0; i<MaxCentBin; i++) {
    TString baseName = Form("EP_%d-%d_%02d",mincent[i],maxcent[i],(Int_t)(10.*emax));
    TString fileName =; 
    if(deta!=deta2) {
      if(mineta2==maxeta2 && maxeta1<0) {
	baseName+="_NegEta";
      } else if(mineta2==maxeta2 && mineta1>0) {
	baseName+="_PosEta";
      } else {
	baseName+="_Asym";
      }
    }
    if(tag.Length()>0) {
      baseName+="-";
      baseName+=tag.Data();
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
    file<<Form("%d-%d     (Npart = %5.1f)",mincent[i],maxcent[i],NpartBin->GetBinContent(i+1))<<endl;
    Double_t * xxx = g[i]->GetX();
    Double_t * yyy = g[i]->GetY();  
    Double_t * yyyerr = g[i]->GetEY();
    for(int j = 0; j< g[i]->GetN();  j++ ) {
      if(yyy[j]<0.0001) continue;
      file<<setprecision(3)<<xxx[j]<<"\t"<<setprecision(3)<<yyy[j]<<"\t"<<yyyerr[j]<<endl;
    }
    file<<tag.Data()<<endl<<endl;
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
    gpt[icent]->Draw("p");
    leg2->AddEntry(g[icent],  Form("%d-%d     (%5.1f)",mincent[icent],maxcent[icent],NpartBin->GetBinContent(icent+1)),"lp");
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
  intV2FileName+=Form("%03d_%03d_",(int)(10.*mineta1),(int)(10.*maxeta1));
  if(maxeta2>mineta2) intV2FileName+=Form("%03d_%03d_",(int)(10.*mineta2),(int)(10.*maxeta2));
  intV2FileName+=tag.Data();
  intV2FileName+=".txt";
  std::ofstream fileInt;
  fileInt.open(intV2FileName.Data());
  double err=0;
  for(int icent = 0; icent<MaxCentBin; icent++) {
    double v2int = IntegralV2(pt[icent], g[icent],gpt[icent],0.3,3.,err);
    cout<<"Integral v2: "<<mincent[icent]<<"-"<<maxcent[icent]<<" = "<<v2int<<" +/- "<<err<<endl;
    fileInt<<setprecision(3)<<aveta1<<"\t"<<setprecision(3)<<v2int<<"\t"<<setprecision(3)<<err<<"\t"<<setprecision(3)<<mincent[icent]<<"\t"<<setprecision(3)<<maxcent[icent]<<"\t"<<endl;
  }
  fileInt.close();
}

double IntegralV2(TH1D * pt,TGraphErrors * g,TGraphErrors * gpt, double ptmin, double ptmax, double &err){
  double ret = 0;
  double norm = 0;
  double norm2 = 0;
  double * gx;
  double * gy;
  double * gey;
  double * gptx;
  double * gpty;

  gx = g->GetX();
  gy = g->GetY();
  gey = g->GetEY();
  gptx = gpt->GetX();
  gpty = gpt->GetY();


  for(Int_t i = 0; i<g->GetN(); i++) {
    int ptbin = pt->FindBin(gx[i]);
    double dpt = pt->GetBinWidth(ptbin);
    if(pt->GetBinLowEdge(ptbin)>=ptmin && pt->GetBinLowEdge(ptbin)<ptmax) {
      norm2+=gpty[i]*dpt;
    }
  }

  err = 0;
  for(Int_t i = 0; i<g->GetN(); i++) {
    int ptbin = pt->FindBin(gx[i]);
    double dpt = pt->GetBinWidth(ptbin);
    if(pt->GetBinLowEdge(ptbin)>=ptmin && pt->GetBinLowEdge(ptbin)<ptmax) {
      ret += gy[i]*gpty[i]*dpt;
      err+=gey[i]*gey[i]*dpt*dpt*gpty[i]*gpty[i]/(norm2*norm2);
      //cout<<i<<" "<<gey[i]<<" "<<dpt<<" "<<gpty[i]<<" "<<norm2<<" "<<gey[i]*dpt*gpty[i]/norm2<<endl;
      norm+=gpty[i]*dpt;
    }
  }
  err=sqrt(err);
  ret/=norm;

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
      if(xx[npt]>minpt) ++npt;
    }
  }
  TGraphErrors * GenV2 = new TGraphErrors(npt,xx,yy,xxerr,yyerr);
  GenV2->SetMarkerStyle(marker);
  GenV2->SetMarkerColor(color);
  GenV2->SetLineColor(color);

  return GenV2;
}
