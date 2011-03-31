#include <iomanip>
#include "TROOT.h"
#include "TSystem.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TVirtualPad.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLine.h"
#include "TList.h"
#include "TString.h"
#include "TVector.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>

Double_t PowerFunc(Double_t * x, Double_t* par) {
  Double_t f;
  Double_t n = par[0];
  Double_t p0 = par[1];
  Double_t dNdy = par[2];
  f = ((n - 1.)*(n - 2.))/(2.*TMath::Pi()*p0*p0);
  
  f = f*dNdy;
  
  f = f*pow((1. + fabs(x[0])/p0),-n);
  return f;
}
Double_t PowerFuncScaled(Double_t * x, Double_t* par) {
  Double_t f;
  Double_t n = par[0];
  Double_t p0 = par[1];
  Double_t dNdy = par[2];
  f = ((n - 1.)*(n - 2.))/(2.*TMath::Pi()*p0*p0);
  
  f = f*dNdy;
  
  f = f*pow((1. + fabs(x[0])/p0),-n);
  f*=2.*TMath::Pi()*x[0];
  return f;
}

Double_t ExpFunc(Double_t * x, Double_t * par) {
  Double_t f;
  Double_t T = par[0];
  Double_t m = par[1];
  Double_t dNdy = par[2];
  
  Double_t mt = TMath::Sqrt( m*m + x[0]*x[0] );
  f = 1./(2.*TMath::Pi()*T*(T+m));
  
  f = f*dNdy;
  
  f = f*TMath::Exp(-(mt-m)/T);
  
  return f;
  
}

Double_t GausFunc(Double_t * x, Double_t * par) {
  Double_t f;
  Double_t sig = par[0];
  Double_t dNdy = par[1];
  
  Double_t A = dNdy/(TMath::TwoPi() * pow(sig,2));
  f = A * TMath::Exp(-x[0]*x[0]/(2*sig*sig));
  
  return f;
  
}

//                   0-5   5-10  10-15    15-20    20-25    25-30     30-35   35-40    40-50  50-60   60-70     70-80 80-90 90-100
int mincent[14] ={     0,     5,    10,      15,      20,      25,      30,     35,       40,    50,     60,       70,   80,    90};
int maxcent[14] ={     5,    10,    15,      20,      25,      30,      35,     40,       50,    60,     70,       80,   90,   100};
int markers[14]= {    22,    21,    23,      24,      25,      30,      20,      3,       29,    28,     26,       27,    3,     5};
int colors[14] = {kBlack, kBlue,  kRed, kCyan+2, kViolet, kSpring, kOrange, kRed+4, kAzure+9, kTeal, kRed-4,kYellow+2,    0,     0};
TFile * tf;

TGraphErrors * GenEtaInt(Int_t minCent=5, Int_t maxCent=10, Double_t detailEtaMin = -0.8, Double_t detailEtaMax = 0.8, TGraphErrors * &v = 0,
			 TGraphErrors *& midSpec, TGraphErrors *& midV2, TGraphErrors *& midV2W, 
			 Double_t &midv2Int=0, Double_t & midv2eInt=0);
TGraphError * GetInt(Double_t minEta, Double_t maxEta, Int_t minCent, Int_t maxCent, Double_t & sInt, Double_t & seInt, Double_t & v2Int, Double_t & v2eInt, TGraphErrors *& sp, TGraphErrors *& v2);
Double_t minX(TF1 * f) {
  Double_t xmin = 0;
  while(f->Eval(xmin)<0) xmin+=0.001;
  return xmin;
}
Double_t XMIN = 0.2;
void GenInt() {
  TString tag = "dzerr14_chi280";
  tf = new TFile(Form("EPSpectra_%s.root",tag.Data()));
  TFile * tfcum = new TFile("Cumul_dNdEta.root");
  Double_t detailEtaMin = -0.8;
  Double_t detailEtaMax = 0.8;
  TGraphErrors * gv2[12];
  TGraphErrors * gSpec[12];
  TGraph * gSpecCum[12];
  TGraphErrors * midV2[12];
  TGraphErrors * midV2W[12];
  TGraphErrors * midSpec[12];
  for(int i = 0; i<12; i++) {
    gSpecCum[i] = (TGraph *) tfcum->Get(Form("dndeta_%d-%d",mincent[i],maxcent[i]));
  }
  TF1 * midFit[12];
  TF1 * midFitR[12];
  Double_t midV2Int[12];
  Double_t midV2EInt[12];
  TLatex * lint[12];
  TLatex * lintFit[12];
  TLatex * lintFull[12];
  TLatex * fullcor[12];
  TLatex * lcent[12];
  for(int i = 0; i<12; i++) {
    gSpec[i]   = GenEtaInt(mincent[i],maxcent[i],detailEtaMin,detailEtaMax,gv2[i],midSpec[i],midV2[i],midV2W[i],midV2Int[i],midV2EInt[i]);
  }
  TCanvas * cv2 = new TCanvas(Form("v2Int_%s",tag.Data()),Form("v2Int_%s",tag.Data()),900,700);
  TH1D * hv2 = new TH1D("hv2","",100,-4,6);
  hv2->SetMaximum(0.12);
  hv2->SetStats(kFALSE);
  hv2->SetXTitle("#eta");
  hv2->SetYTitle("v_{2}[EP]");
  hv2->Draw();
  gPad->SetGrid(1,1);
  TLegend * leg = new TLegend(0.7,0.4,0.88,0.88,"Centrality");
  for(int i = 0; i< 12; i++) {
    gv2[i]->SetMarkerStyle(markers[i]);
    gv2[i]->SetMarkerColor(colors[i]);
    gv2[i]->SetLineColor(colors[i]);
    gv2[i]->Draw("p");
    leg->AddEntry(gv2[i],Form("%d - %d%c",mincent[i],maxcent[i],'%'),"p");
  }
  leg->Draw();
  TCanvas * cspec = new TCanvas(Form("cspec_%s",tag.Data()),Form("cspec_%s",tag.Data()),900,700);
  TH1D * hspec = new TH1D("hspec","",100,-4,6);
  hspec->SetMaximum(1200);
  hspec->SetStats(kFALSE);
  hspec->SetXTitle("#eta");
  hspec->SetYTitle("dN/d#eta");
  hspec->Draw();
  gPad->SetGrid(1,1);
  TLegend * leg2 = new TLegend(0.7,0.4,0.88,0.88,"Centrality");
  for(int i = 0; i< 12; i++) {
    gSpec[i]->SetMarkerStyle(markers[i]);
    gSpec[i]->SetMarkerColor(kBlue);
    gSpec[i]->SetLineColor(colors[i]);
    gSpec[i]->Draw("p");
    leg2->AddEntry(gSpec[i],Form("%d - %d%c",mincent[i],maxcent[i],'%'),"p");

    if(gSpecCum[i]) {
      gSpecCum[i]->SetMarkerStyle(markers[i]);
      gSpecCum[i]->SetMarkerColor(kRed);
      gSpecCum[i]->Draw("p");
    }
  }
  leg2->Draw();
  TLatex * EP = new TLatex(3.2,280,"EP");
  EP->SetTextColor(kBlue);
  EP->Draw();
  TLatex * Cum = new TLatex(3.2,210,"Cumulant");
  Cum->SetTextColor(kRed);
  Cum->Draw();
  cv2->Print(Form("v2Int_%s.png",tag.Data()),"png");
  cspec->Print(Form("dNdEta_%s.png",tag.Data()),"png");
  cv2->Print(Form("v2Int_%s.pdf",tag.Data()),"pdf");
  cspec->Print(Form("dNdEta_%s.pdf",tag.Data()),"pdf");

  TString detailName = Form("detailV2Int_%s_%04.1f_%04.1f_ptMin_%03.1lf",tag.Data(),detailEtaMin,detailEtaMax,XMIN);
  TCanvas * cv2w = new TCanvas(detailName.Data(),detailName.Data(),1200,850);
  cv2w->Divide(4,3);
  TGraphErrors * midSpecNorm[12];
  TGraphErrors * invSpec[12];
  TH1D * hmidSpecNorm[12];
  Double_t yy[50];
  Double_t eyy[50];
  Double_t yyinv[50];
  Double_t eyyinv[50];
  Double_t dndy[12]={1601,1294,1070,870,720,590,470,380,261,149,76,35};
  TF1 *sfit[12];
  TF1 * sfitScaled[12];
  TF1 * sfitScaledFull[12];
  for(int ic = 0; ic<12; ic++) {
    cv2w->cd(ic+1);
    Double_t * xSpec = midSpec[ic]->GetX();
    Double_t * ySpec = midSpec[ic]->GetY();
    Double_t * eySpec = midSpec[ic]->GetEY();
    Double_t * xV2W = midV2W[ic]->GetX();
    Double_t * yV2W = midV2W[ic]->GetY();
    Double_t * eyV2W = midV2W[ic]->GetEY();
    Double_t * xV2 = midV2[ic]->GetX();
    Double_t * yV2 = midV2[ic]->GetY();
    Double_t * eyV2 = midV2[ic]->GetEY();
    Double_t sum = 0;
    for(int i = 0; i< midSpec[ic]->GetN(); i++) sum+=ySpec[i];
    for(int i = 0; i< midSpec[ic]->GetN(); i++) { 
      yy[i]=ySpec[i]/sum;
      eyy[i]=eySpec[i]/sum;
      yyinv[i]=ySpec[i]/xSpec[i]/2./TMath::Pi();
      eyyinv[i]=eySpec[i]/xSpec[i]/2./TMath::Pi();
      eyyinv[i]+=TMath::Sqrt(eyyinv[i]*eyyinv[i]+pow(0.05*yyinv[i],2));
      cout<<xSpec[i]<<" "<<ySpec[i]<<" "<<eySpec[i]<<" "<<yyinv[i]<<" "<<eyyinv[i]<<endl;
    }
    sfit[ic] = new TF1("sfit",PowerFunc,0,3,3);
    sfit[ic]->SetParNames("n","p0","dN/dy");
    sfit[ic]->SetParameters(28,5,dndy[ic]);
    sfit[ic]->FixParameter(2,dndy[ic]);
    sfitScaled[ic] = new TF1(Form("sfitScaled_%d",ic),PowerFuncScaled,0.3,3.,3);
    sfitScaledFull[ic] = new TF1(Form("sfitScaledFull_%d",ic), PowerFuncScaled,0.,3.,3);
    invSpec[ic] = new TGraphErrors(midSpec[ic]->GetN(),xSpec,yyinv,0,eyyinv);
    invSpec[ic]->Fit(sfit[ic],"RN");
    sfitScaled[ic]->SetParameters( sfit[ic]->GetParameter(0),sfit[ic]->GetParameter(1),sfit[ic]->GetParameter(2)/sum);
    sfitScaledFull[ic]->SetParameters( sfit[ic]->GetParameter(0),sfit[ic]->GetParameter(1),sfit[ic]->GetParameter(2)/sum);
    cout<<"Chisq: "<<sfit[ic]->GetChisquare()<<"    Ndf: "<<sfit[ic]->GetNDF()<<endl;
    midSpecNorm[ic] = new TGraphErrors(midSpec[ic]->GetN(),xSpec,yy,0,eyy);
    TString hspecname = Form("hSpec%d_%d",mincent[ic],maxcent[ic]);
    TString hspectitle=Form("%d - %d\%",mincent[ic],maxcent[ic]);
    hmidSpecNorm[ic]=new TH1D(hspecname.Data(),hspectitle.Data(),100,0,4);
    hmidSpecNorm[ic]->SetMaximum(0.4);
    hmidSpecNorm[ic]->SetMinimum(0.0);
    hmidSpecNorm[ic]->SetXTitle("p_{T} (GeV/c)");
    hmidSpecNorm[ic]->SetStats(kFALSE);
    hmidSpecNorm[ic]->Draw();
    midSpecNorm[ic]->SetMarkerStyle(markers[ic]);
    midSpecNorm[ic]->SetMarkerColor(colors[ic]);
    midSpecNorm[ic]->SetLineStyle(colors[ic]);
    midSpecNorm[ic]->Draw("p");
    sfitScaled[ic]->SetLineColor(colors[ic]);
    sfitScaled[ic]->Draw("same");
    sfitScaledFull[ic]->SetLineColor(colors[ic]);
    sfitScaledFull[ic]->SetLineWidth(0.8);
    sfitScaledFull[ic]->Draw("same");
    for(Double_t x = 0.01; x<3; x+=0.1) cout<<x<<" "<<sfitScaled[ic]->Eval(x)<<endl;
    midV2[ic]->Draw("lp");
    midV2W[ic]->Draw("p");
    midFit[ic] = new TF1(Form("fit_%d",ic),"pol6",0.0,3.2);
    midFitR[ic] = new TF1(Form("fitR_%d",ic),"pol6",0.3,3.0);
    midFit[ic]->SetFillStyle(2001);
    midFit[ic]->SetFillColor(colors[ic]);
    midFit[ic]->SetLineWidth(1);
    midFit[ic]->FixParameter(0,0.);
    midV2W[ic]->Fit(midFit[ic],"R");
    midV2W[ic]->GetFunction(Form("fit_%d",ic))->SetFillStyle(3004);
    midV2W[ic]->GetFunction(Form("fit_%d",ic))->SetFillColor(colors[ic]);
    midV2W[ic]->GetFunction(Form("fit_%d",ic))->Draw("SAME FC");
    midFitR[ic]->SetParameters(midV2W[ic]->GetFunction(Form("fit_%d",ic))->GetParameters());
    midFitR[ic]->SetFillStyle(3002);
    midFitR[ic]->SetFillColor(colors[ic]);
    midFitR[ic]->Draw("SAME FC");
    TF1 * tmp = midV2W[ic]->GetFunction(Form("fit_%d",ic));
    Double_t xmin=XMIN;
    if(minX(tmp)>xmin)  xmin = minX(tmp);
    cout<<"xmin: "<<minX(tmp)<<endl;
    Double_t intRange = tmp->Integral(0.3,3.0)/2.;
    Double_t intFull =  tmp->Integral(xmin,3.0)/2.;
    intFull *= sfitScaledFull[ic]->Integral(0.3,3.0)/sfitScaledFull[ic]->Integral(XMIN,3.);
    cout<<"int: "<<intRange<<" "<<intFull<<" "<<midV2Int[ic]<<" "<<midV2EInt[ic]<<endl;
    lint[ic] = new TLatex(0.8,0.32,Form(    "v_{2}[EP;0.3-3.0 GeV/c] = %6.3f",midV2Int[ic]));
    lint[ic]->Draw();
    lintFit[ic] = new TLatex(0.8,0.29,Form( "v_{2}[EP; fit in range] = %6.3f",intRange));
    lintFit[ic]->Draw();
    lintFull[ic] = new TLatex(0.8,0.26,Form("v_{2}[EP;%3.1lf-3.0 GeV/c] = %6.3f",XMIN,intFull));
    lintFull[ic]->Draw();
    fullcor[ic] = new TLatex(0.8,0.23,Form("Missing: %6.1f%c",100.*(intRange-intFull)/intFull,'%'));
    fullcor[ic]->Draw();
    lcent[ic] = new TLatex(0.3,0.36,Form("%d - %d%c",mincent[ic],maxcent[ic],'%'));
    lcent[ic]->Draw();
  }
  cv2w->Print(Form("%s.png",detailName.Data()),"png");
  cv2w->Print(Form("%s.pdf",detailName.Data()),"pdf");
}

TGraphErrors * GenEtaInt(Int_t minCent, Int_t maxCent, Double_t detailEtaMin, Double_t detailEtaMax,TGraphErrors *&v2, TGraphErrors *& midSpec, TGraphErrors *&midV2,
			 TGraphErrors *& midV2W, Double_t &midv2Int, Double_t & midv2eInt) {
  Double_t etabins[]={-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0.,0.4,0.8,1.2,1.6,2.0,2.4};
  Double_t yint[12];
  Double_t yeint[12];
  Double_t yv2[12];
  Double_t yev2[12];
  Double_t xint[12];
  TGraphErrors * etaSpec[12];
  TGraphErrors * etaV2[12];
  TGraphErrors * etaV2W[12];
  cout<<minCent<<" <= cent < "<<maxCent<<endl;
  for(Int_t ie = 0; ie<12; ie++) {
    Double_t sInt=0;
    Double_t v2Int=0;
    Double_t seInt = 0;
    Double_t v2eInt = 0;
    etaV2W[ie]=GetInt(etabins[ie],etabins[ie+1],minCent,maxCent,sInt,seInt, v2Int, v2eInt, etaSpec[ie], etaV2[ie]);
    xint[ie] = etabins[ie]+0.2;
    yint[ie] = sInt;
    yeint[ie] = seInt;
    yv2[ie] = v2Int;
    yev2[ie]=v2eInt;
  }
  sInt = 0;
  midv2Int = 0;
  seInt = 0;
  midv2eInt = 0;
  midV2W = GetInt(detailEtaMin,detailEtaMax,minCent,maxCent,sInt,seInt,midv2Int,midv2eInt,midSpec,midV2);

  cout<<"====="<<endl;
  TGraphErrors * ret = new TGraphErrors(12,xint,yint,0,yeint);
  v2 = new TGraphErrors(12,xint,yv2,0,yev2);
  return ret;
}

TGraphErrors *  GetInt(Double_t minEta, Double_t maxEta, Int_t minCent, Int_t maxCent,Double_t &sInt, Double_t & seInt, Double_t &v2Int, Double_t &v2eInt, TGraphErrors *& sp, TGraphErrors *& v2) {
  TH1D * ptbins = (TH1D *) tf->Get("ptbins");
  TString sname = Form("%04.1f_%04.1f/dNdPt_%d_%d",minEta,maxEta,minCent,maxCent);
  if(maxEta-minEta>1.) sname = Form("-0.8_00.0_00.0_00.8/dNdPt_%d_%d",minCent,maxCent);
  TString v2name = Form("%04.1f_%04.1f/v2_%d_%d",minEta,maxEta,minCent,maxCent);
  if(maxEta-minEta>1.) v2name = Form("-0.8_00.0_00.0_00.8/v2_%d_%d",minCent,maxCent);
  Bool_t pr = kFALSE;
  if(maxEta-minEta>1.) pr = kFALSE;
  cout<<sname.Data()<<endl;
  cout<<v2name.Data()<<endl;
  TGraphErrors * gs = tf->Get(sname.Data());
  TGraphErrors * gv2 = tf->Get(v2name.Data());
  sp = gs;
  v2 = gv2;
  Double_t * x = gs->GetX();
  Double_t * y = gs->GetY();
  Double_t * ey = gs->GetEY();
  Double_t * yv2 = gv2->GetY();
  Double_t * yv2e = gv2->GetEY();
  sInt=0;
  seInt = 0;
  v2Int = 0;
  v2eInt = 0;
  Double_t xx[40];
  Double_t yy[40];
  Double_t eyy[40];
  Int_t n = 0;
  if(pr)cout<<"GetInt for cent: "<<minCent<<" "<<maxCent<<"  and eta range: "<<minEta<<" "<<maxEta<<endl;
  for(Int_t i = 0; i<gs->GetN(); i++) {
    Double_t pt = x[i];
    Int_t ptbin = ptbins->FindBin(pt);
    Double_t wpt = ptbins->GetBinWidth(ptbin);
    if(pt>=0.3&&pt<3.0) {
      sInt+=y[i]*wpt;
      seInt+=pow(ey[i]*wpt,2);
      v2Int+=yv2[i]*y[i]*wpt;
      v2eInt+=pow(yv2e[i]*y[i]*wpt,2);
      xx[n]=pt;
      yy[n] = yv2[i]*y[i];
      eyy[n++] = yv2e[i]*y[i]; 
      if(pr)cout<<x[i]<<" "<<y[i]<<"+/-"<<ey[i]<<" -----  "<<yv2[i]<<" "<<yv2e[i]<<endl;
      if(pr)cout<<"=="<<xx[n-1]<<" "<<yy[n-1]<<"+/-"<<eyy[n-1]<<"   wpt: "<<wpt<<endl;
   }
  }
  seInt=sqrt(seInt);
  v2eInt=sqrt(v2eInt);
  v2Int/=sInt;
  v2eInt/=sInt;
  if(pr) cout<<"sInt: "<<sInt<<"   v2Int: "<<v2Int<<endl;
  for(int i = 0; i<n; i++) {
    yy[i]/=sInt/2.;
    eyy[i]/=sInt/2.;
  }
  xx[n] = 3.2;
  yy[n] = 0.;
  eyy[n] = 0.001;
  cout<<setprecision(4)<<minEta<<" "<<maxEta<<"   dNdEta: "<<sInt<<"+/-"<<seInt<<"    v2: "<<v2Int<<" "<<v2eInt<<endl;
  TGraphErrors * gw = new TGraphErrors(n+1,xx,yy,0,eyy);
  gw->SetMarkerStyle(30);
  gw->SetMarkerColor(kBlue);
  gw->SetLineColor(kBlue);
  return gw;
}
