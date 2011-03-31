#include <stdlib.h>
#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TString.h"
#include "TCanvasImp.h"
#include "TGraphErrors.h"
#include <stdio.h>
#include <iostream>
TGraphErrors * InputType1(TString fileName, Int_t & Num);
TGraphErrors * InputType1(int  minCent, int maxCent, double eta, TString tag, TString EP1, TString EP2, Int_t & Num) {
  TString name = Form("results/v2pt_EP_%d-%d_%02d-_%s_%s_%s.txt",minCent,maxCent,(int)(10*eta),tag.Data(),EP1.Data(),EP2.Data());
  if(EP2.Length()==0) name = Form("results/v2pt_EP_%d-%d_%02d_%s_%s.txt",minCent,maxCent,(int)(10*eta+0.1),tag.Data(),EP1.Data());
  return InputType1(name,Num);
}
TGraphErrors * spec; 
static const double ptbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0};

void CompareEta() {

  Int_t minCent = 10;
  Int_t maxCent = 15;
  int ptbin = 7;
  TGraphErrors * gV2PosEP[20];
  TGraphErrors * gSpecPosEP[20];
  double etaPos[20];
  double etaNeg[20];
  TString EP1 = "etHF";
  int npos = 0;
  Int_t nPts1= 0;
  for(double eta = -2.4; eta<=1.0; eta = eta+0.2) {
    etaPos[npos] = eta+0.1;
    TString tag = "NegEta-_hiGoodMergedTracks_Flow_Skim_Run2010-v5_dz5Flat_-5to5";
    if(eta>=-0.01) tag = "PosEta-_hiGoodMergedTracks_Flow_Skim_Run2010-v5_dz5Flat_-5to5";
    double etaval = eta;
    if(eta>=-0.0001) etaval = eta+0.2;
    gV2PosEP[npos] = InputType1(minCent,maxCent,fabs(etaval),tag,EP1,"",nPts1);
    gSpecPosEP[npos] = spec;
    ++npos;
  }
  TGraphErrors * gV2NegEP[20];
  TGraphErrors * gSpecNegEP[20];
  TString EP1 = "etHF";
  int nneg = 0;
  nPts1= 0;
  for(double eta = -1.2; eta<=2.2; eta = eta+0.2) {
    etaNeg[nneg] = eta+0.1;
    TString tag = "NegEta-_hiGoodMergedTracks_Flow_Skim_Run2010-v5_dz5Flat_-5to5";
    if(eta>=-0.01) tag = "PosEta-_hiGoodMergedTracks_Flow_Skim_Run2010-v5_dz5Flat_-5to5";
    double etaval = eta;
    if(eta>=-0.0001) etaval = eta+0.2;
    gV2NegEP[nneg] = InputType1(minCent,maxCent,fabs(etaval),tag,EP1,"",nPts1);
    gSpecNegEP[nneg] = spec;
    ++nneg;
  }

  TCanvas * c1 = new TCanvas("c1","c1",1000,600);
  c1->Divide(2);
  c1->cd(1);
  gPad->SetGrid(1,1);
  double xpos[20];
  double ypos[20];
  double dypos[20];
  double yypos[20];
  double dyypos[20];
  double xneg[20];
  double yneg[20];
  double dyneg[20];
  double yyneg[20];
  double dyyneg[20];
  for(int i = 0; i<npos; i++) {
    Double_t x;
    Double_t y;
    double dy;
    gV2PosEP[i]->GetPoint(ptbin,x,y);
    dy = gV2PosEP[i]->GetErrorY(ptbin);
    xpos[i] = etaPos[i];
    ypos[i] = y;
    dypos[i] = dy;
    gSpecPosEP[i]->GetPoint(ptbin,x,y);
    yypos[i] = y;
    dyypos[i] = gSpecPosEP[i]->GetErrorY(ptbin);

  }
  for(int i = 0; i<nneg; i++) {
    Double_t x;
    Double_t y;
    double dy;
    xneg[i] = gV2NegEP[i]->GetPoint(ptbin,x,y);
    dy = gV2NegEP[i]->GetErrorY(ptbin);
    xneg[i] = etaNeg[i];
    yneg[i] = y;
    dyneg[i] = dy;
    gSpecNegEP[i]->GetPoint(ptbin,x,y);
    yyneg[i] = y;
    dyyneg[i] = gSpecNegEP[i]->GetErrorY(ptbin);
  }
 
  double ref = ( yypos[10]+yypos[11]+ yypos[12]+yyneg[5]+yyneg[6]+yyneg[7])/6.;
  double v2ref = (ypos[10]+ypos[11]+  ypos[12]+ yneg[5]+ yneg[6]+yneg[7])/6.;
  for(int i = 0; i<npos; i++) {
    ypos[i]/=v2ref;
    dypos[i]/=v2ref;
    yypos[i]/=ref;
    dyypos[i]/=ref;
  }
  for(int i = 0; i<nneg; i++) {
    yneg[i]/=v2ref;
    dyneg[i]/=v2ref;
    yyneg[i]/=ref;
    dyyneg[i]/=ref;
  }
  TGraphErrors * ptp = new TGraphErrors(npos,xpos,ypos,0,dypos);
  TGraphErrors * ptsp = new TGraphErrors(npos,xpos,yypos,0,dyypos);
  ptp->SetMarkerStyle(21);
  ptp->SetMarkerColor(kRed);
  ptsp->SetMarkerStyle(21);
  ptsp->SetMarkerColor(kRed);
  TH1D * href = new TH1D("v2Eta",Form("v_{2} vs. #eta (p_{T} = %4.2lf GeV/c)",(ptbins[ptbin]+ptbins[ptbin+1])/2.),200,-3,3);
  href->SetMaximum(1.4);
  href->SetMinimum(0.6);
  href->SetStats(kFALSE);
  href->SetXTitle("#eta");
  href->SetYTitle("Ratio to mid-#eta");
  href->Draw();
  ptp->Draw("p");
  TGraphErrors * ptn = new TGraphErrors(nneg,xneg,yneg,0,dyneg);
  TGraphErrors * ptsn = new TGraphErrors(nneg,xneg,yyneg,0,dyyneg);
  ptn->SetMarkerStyle(25);
  ptn->SetMarkerColor(kBlue);
  ptsn->SetMarkerStyle(25);
  ptsn->SetMarkerColor(kBlue);
  ptn->Draw("p");
  TLatex * lab = new TLatex(0.0,0.65,Form("Centrality: %d - %d",minCent,maxCent));
  lab->Draw();
  TLegend * leg = new TLegend(0.55,0.7,0.88,0.88);
  leg->AddEntry(ptp,"etHF","p");
  leg->AddEntry(ptn,"etHF","p");
  leg->Draw();
  TLine * line = new TLine(-2.6,1.,2.6,1.);
  line->SetLineColor(kGreen+1);
  line->SetLineStyle(2);
  line->SetLineWidth(2.0);
  line->Draw();
  c1->cd(2);
  TH1D * href2 = new TH1D("SpecEta","(1/N_{evt})d^{2}N/dp_{T}d#eta/(#eta=0 ref) vs. #eta",200,-3,3);
  href2->SetMaximum(1.4);
  href2->SetMinimum(0.6);
  gPad->SetGrid(1,1);
  href2->SetStats(kFALSE);
  href2->SetXTitle("#eta");
  href2->SetYTitle("Ratio to mid-#eta");
  href2->Draw();

  ptsp->Draw("p");
  ptsn->Draw("p");
  line->Draw();
  c1->Print(Form("EtaDists_%d_%d_%d.png",minCent,maxCent,(Int_t)((ptbins[ptbin]+ptbins[ptbin+1])*50.)),"png");
}

TGraphErrors * InputType1(TString fileName,  Int_t &Num){
  FILE * fin = fopen(fileName.Data(),"r");
  Double_t x[100];
  Double_t y[100];
  Double_t dy[100];
  Double_t dNdPt[100];
  Double_t dNdPtErr[100];
  Num = 1;
  if(fin==NULL) {
    cout<<"File not found: "<<fileName.Data()<<endl;
    return 0;
  }
  char buf[80];
  fgets(buf,80,fin);
  while(fgets(buf,80,fin)!=NULL) {
    if(buf[0]=='_') continue;
    if(sscanf(buf,"%lf\t%lf\t%lf\t%lf\t%lf",&x[Num],&y[Num],&dy[Num],&dNdPt[Num],&dNdPtErr[Num]) == 5) {
      ++Num;
    }; 

  }
  if(Num==0) return 0;
  TGraphErrors * g = new TGraphErrors(Num,x,y,0,dy);
  spec = new TGraphErrors(Num,x,dNdPt,0,dNdPtErr);
  g->SetName(fileName.Data());
  g->SetTitle(fileName.Data());
  return g;
}
