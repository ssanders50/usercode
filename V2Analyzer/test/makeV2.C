#include "/net/hisrv0001/home/sanders/CMSSW_3_9_1/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include <iomanip>
#include <fstream>
using namespace std;
TFile * tf;
void makeV2(TString tag = "Hydjet_Bass"){
  tf = new TFile("/net/hisrv0001/home/sanders/CMSSW_3_9_1/src/V2Analyzer/V2Analyzer/data/rpflat_combined.root");
  TCanvas * can = new TCanvas(tag.Data(),tag.Data(),800,500);
  // Step 1.   Generate EP Resolutions for pos and neg track EP
  Double_t centBins[]={0,5,10,20,30,40,50,60,70,80,90,100};
  TH1D * resNeg = new TH1D("resNeg","resNeg",11,centBins);
  TH1D * resPos = new TH1D("resPos","resPos",11,centBins);
  TH1D * hNpartBin = new TH1D("hNpartBin","hNpartBin",11,centBins);
  TH2D * EPCorr[20];
  TH2D * hpt = (TH2D *) tf->Get("v2analyzer/pt")->Clone("hpt");
  TH2D * dNdPt = (TH2D *) tf->Get("v2analyzer/ptCnt")->Clone("dNdPt");
  hpt->Divide(dNdPt);
  for(int i = 1; i<= dNdPt->GetNbinsX(); i++ ) 
    for(int j = 1; j<=20; j++) {
      dNdPt->SetBinContent(i,j, dNdPt->GetBinContent(i,j)/
			   dNdPt->GetXaxis()->GetBinWidth(i)/2.);
    }
  TH1D * NpartBin = (TH1D *) tf->Get("v2analyzer/NpartBin");
  TH1D * NpartBinCnt = (TH1D *) tf->Get("v2analyzer/NpartBinCnt");

  EPCorr[0] = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_0")->Clone("EPCorr_0_5");
  EPCorr[0]->Divide((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_0"));
  hNpartBin->SetBinContent(1,NpartBin->GetBinContent(1)/NpartBinCnt->GetBinContent(1));

  EPCorr[1] = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_1")->Clone("EPCorr_5_10");
  EPCorr[1]->Divide((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_1"));
  hNpartBin->SetBinContent(2,NpartBin->GetBinContent(2)/NpartBinCnt->GetBinContent(2));

  EPCorr[2] =    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_2")->Clone("EPCorr_10_20");
  EPCorr[2]->Add((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_3"));
  TH2D * tmp = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_2");
  tmp->Add(    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_3"));
  EPCorr[2]->Divide(tmp);
  hNpartBin->SetBinContent(3,(NpartBin->GetBinContent(3)+NpartBin->GetBinContent(4))/
			   (NpartBinCnt->GetBinContent(3)+NpartBinCnt->GetBinContent(4)));

  EPCorr[3] =    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_4")->Clone("EPCorr_20_30");
  EPCorr[3]->Add((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_5"));
  tmp = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_4");
  tmp->Add(    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_5"));
  EPCorr[3]->Divide(tmp);
  hNpartBin->SetBinContent(4,(NpartBin->GetBinContent(5)+NpartBin->GetBinContent(6))/
			   (NpartBinCnt->GetBinContent(5)+NpartBinCnt->GetBinContent(6)));

  EPCorr[4] =    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_6")->Clone("EPCorr_30_40");
  EPCorr[4]->Add((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_7"));
  tmp = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_6");
  tmp->Add(    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_7"));
  EPCorr[4]->Divide(tmp);
  hNpartBin->SetBinContent(5,(NpartBin->GetBinContent(7)+NpartBin->GetBinContent(8))/
			   (NpartBinCnt->GetBinContent(7)+NpartBinCnt->GetBinContent(8)));


  EPCorr[5] =    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_8")->Clone("EPCorr_40_50");
  EPCorr[5]->Add((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_9"));
  tmp = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_8");
  tmp->Add(    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_9"));
  EPCorr[5]->Divide(tmp);
  hNpartBin->SetBinContent(6,(NpartBin->GetBinContent(9)+NpartBin->GetBinContent(10))/
			   (NpartBinCnt->GetBinContent(9)+NpartBinCnt->GetBinContent(10)));

  EPCorr[6] =    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_10")->Clone("EPCorr_50_60");
  EPCorr[6]->Add((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_11"));
  tmp = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_10");
  tmp->Add(    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_11"));
  EPCorr[6]->Divide(tmp);
  hNpartBin->SetBinContent(7,(NpartBin->GetBinContent(11)+NpartBin->GetBinContent(12))/
			   (NpartBinCnt->GetBinContent(11)+NpartBinCnt->GetBinContent(12)));

  EPCorr[7] =    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_12")->Clone("EPCorr_60_70");
  EPCorr[7]->Add((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_13"));
  tmp = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_12");
  tmp->Add(    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_13"));
  EPCorr[7]->Divide(tmp);
  hNpartBin->SetBinContent(8,(NpartBin->GetBinContent(13)+NpartBin->GetBinContent(14))/
			   (NpartBinCnt->GetBinContent(13)+NpartBinCnt->GetBinContent(14)));

  EPCorr[8] =    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_14")->Clone("EPCorr_70_80");
  EPCorr[8]->Add((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_15"));
  tmp = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_14");
  tmp->Add(    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_15"));
  EPCorr[8]->Divide(tmp);
  hNpartBin->SetBinContent(9,(NpartBin->GetBinContent(15)+NpartBin->GetBinContent(16))/
			   (NpartBinCnt->GetBinContent(15)+NpartBinCnt->GetBinContent(16)));

  EPCorr[9] =    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_16")->Clone("EPCorr_80_90");
  EPCorr[9]->Add((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_17"));
  tmp = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_16");
  tmp->Add(    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_17"));
  EPCorr[9]->Divide(tmp);
  hNpartBin->SetBinContent(10,(NpartBin->GetBinContent(17)+NpartBin->GetBinContent(18))/
			   (NpartBinCnt->GetBinContent(17)+NpartBinCnt->GetBinContent(18)));


  EPCorr[10] =    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_18")->Clone("EPCorr_90_100");
  EPCorr[10]->Add((TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelation_19"));
  tmp = (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_18");
  tmp->Add(    (TH2D *) tf->Get("v2analyzer/EPCorrelations/EPCorrelationCnt_19"));
  EPCorr[10]->Divide(tmp);
  hNpartBin->SetBinContent(11,NpartBin->GetBinContent(19)+NpartBin->GetBinContent(20)/
			   (NpartBinCnt->GetBinContent(19)+NpartBinCnt->GetBinContent(20)));

  for(int i = 0; i< 11; i++) {
    Double_t argP = EPCorr[i]->GetBinContent(EvtPTracksPosEtaGap+1,EvtPlaneFromTracksMidEta+1) *
      EPCorr[i]->GetBinContent(EvtPTracksPosEtaGap+1,EvtPTracksNegEtaGap+1)/
      EPCorr[i]->GetBinContent(EvtPlaneFromTracksMidEta+1,EvtPTracksNegEtaGap+1); 
    Double_t argN = EPCorr[i]->GetBinContent(EvtPTracksNegEtaGap+1,EvtPlaneFromTracksMidEta+1) *
      EPCorr[i]->GetBinContent(EvtPTracksPosEtaGap+1,EvtPTracksNegEtaGap+1)/
      EPCorr[i]->GetBinContent(EvtPlaneFromTracksMidEta+1,EvtPTracksPosEtaGap+1); 
    cout<<i<<" "<<argP<<" "<<argN<<endl;
    if(argP<=0 || argN <=0 ) continue;
    resPos->SetBinContent(i+1,sqrt(argP));
    resNeg->SetBinContent(i+1,sqrt(argN));
  }
  TH2D * pos = (TH2D *) tf->Get("v2analyzer/v2/v2Reco/EvtPTracksPosEtaGap/cos_EvtPTracksPosEtaGap_1");
  pos->Divide( (TH2D *) tf->Get("v2analyzer/v2/v2Reco/EvtPTracksPosEtaGap/cnt_EvtPTracksPosEtaGap_1"));
  TH2D * neg = (TH2D *) tf->Get("v2analyzer/v2/v2Reco/EvtPTracksNegEtaGap/cos_EvtPTracksNegEtaGap_2");
  neg->Divide( (TH2D *) tf->Get("v2analyzer/v2/v2Reco/EvtPTracksNegEtaGap/cnt_EvtPTracksNegEtaGap_2"));
  TH1D * pos_0_5 = (TH1D *) pos->ProjectionY("pos_0_5",1,1);
  pos_0_5->Scale(1./resPos->GetBinContent(1));
  TH1D * neg_0_5 = (TH1D *) neg->ProjectionY("neg_0_5",1,1);
  neg_0_5->Scale(1./resNeg->GetBinContent(1));
  TH1D * v2_0_5 = (TH1D *) pos_0_5->Clone("v2_0_5");
  v2_0_5->Add(neg_0_5);
  v2_0_5->Scale(0.5);
  v2_0_5->SetMinimum(0.0);
  v2_0_5->SetMaximum(0.25);
  int marker_0_5 = 22;
  int color_0_5 = kBlack;
  v2_0_5->SetMarkerStyle(marker_0_5);
  neg->SetMinimum(0);
  neg->SetMaximum(0.4);
  Int_t numpt = v2_0_5->GetNbinsX();
  Double_t xx[25];
  Double_t yy[25];
  Double_t xxerr[25];
  Double_t yyerr[25];
  for(int i = 0; i< numpt; i++ ) {
    yy[i] = v2_0_5->GetBinContent(i+1);
    yyerr[i] = v2_0_5->GetBinError(i+1);
    xx[i] = hpt->GetBinContent(i+1,1);
    xxerr[i]=0;
  }
  TGraphErrors * g_0_5 = new TGraphErrors(numpt,xx,yy,xxerr,yyerr);
  g_0_5->SetMarkerStyle(marker_0_5);
  g_0_5->SetMarkerColor(color_0_5);
  g_0_5->SetLineColor(color_0_5);

  Int_t marker_5_10 = 21;
  Int_t color_5_10 = kBlue;
  TH1D * pos_5_10 = (TH1D *) pos->ProjectionY("pos_5_10",2,2);
  pos_5_10->Scale(1./resPos->GetBinContent(2));
  TH1D * neg_5_10 = (TH1D *) neg->ProjectionY("neg_5_10",2,2);
  neg_5_10->Scale(1./resNeg->GetBinContent(2));
  TH1D * v2_5_10 = (TH1D *) pos_5_10->Clone("v2_5_10");
  v2_5_10->Add(neg_5_10);
  v2_5_10->Scale(0.5);
  v2_5_10->SetMinimum(-0.1);
  v2_5_10->SetMaximum(0.4);
  v2_5_10->SetMarkerStyle(marker_5_10);
  v2_5_10->SetMarkerColor(color_5_10);
  v2_5_10->SetLineColor(color_5_10);
  for(int i = 0; i< numpt; i++ ) {
    yy[i] = v2_5_10->GetBinContent(i+1);
    yyerr[i] = v2_5_10->GetBinError(i+1);
    xx[i] = hpt->GetBinContent(i+1,2);
    xxerr[i]=0;
  }
  TGraphErrors * g_5_10 = new TGraphErrors(numpt,xx,yy,xxerr,yyerr);
  g_5_10->SetMarkerStyle(marker_5_10);
  g_5_10->SetMarkerColor(color_5_10);
  g_5_10->SetLineColor(color_5_10);

  int marker_10_20 = 23;
  int color_10_20 = kRed;
  TH1D * pos_10_20 = (TH1D *) pos->ProjectionY("pos_10_20",3,4);
  pos_10_20->Scale(1./resPos->GetBinContent(3));
  TH1D * neg_10_20 = (TH1D *) neg->ProjectionY("neg_10_20",3,4);
  neg_10_20->Scale(1./resNeg->GetBinContent(3));
  TH1D * v2_10_20 = (TH1D *) pos_10_20->Clone("v2_10_20");
  v2_10_20->Add(neg_10_20);
  v2_10_20->Scale(0.25);
  v2_10_20->SetMinimum(-0.1);
  v2_10_20->SetMaximum(0.4);
  v2_10_20->SetMarkerStyle(marker_10_20);
  v2_10_20->SetMarkerColor(color_10_20);
  v2_10_20->SetLineColor(color_10_20);
  for(int i = 0; i< numpt; i++ ) {
    yy[i] = v2_10_20->GetBinContent(i+1);
    yyerr[i] = v2_10_20->GetBinError(i+1);
    xx[i] = (hpt->GetBinContent(i+1,3)+hpt->GetBinContent(i+1,4))/2.;
    xxerr[i]=0;
  }
  TGraphErrors * g_10_20 = new TGraphErrors(numpt,xx,yy,xxerr,yyerr);
  g_10_20->SetMarkerStyle(marker_10_20);
  g_10_20->SetMarkerColor(color_10_20);
  g_10_20->SetLineColor(color_10_20);

  int marker_20_30 = 24;
  int color_20_30 = kCyan+2;
  TH1D * pos_20_30 = (TH1D *) pos->ProjectionY("pos_20_30",5,6);
  pos_20_30->Scale(1./resPos->GetBinContent(4));
  TH1D * neg_20_30 = (TH1D *) neg->ProjectionY("neg_20_30",5,6);
  neg_20_30->Scale(1./resNeg->GetBinContent(4));
  TH1D * v2_20_30 = (TH1D *) pos_20_30->Clone("v2_20_30");
  v2_20_30->Add(neg_20_30);
  v2_20_30->Scale(0.25);
  v2_20_30->SetMinimum(-0.1);
  v2_20_30->SetMaximum(0.4);
  v2_20_30->SetMarkerStyle(marker_20_30);
  v2_20_30->SetMarkerColor(color_20_30);
  v2_20_30->SetLineColor(color_20_30);
  for(int i = 0; i< numpt; i++ ) {
    yy[i] = v2_20_30->GetBinContent(i+1);
    yyerr[i] = v2_20_30->GetBinError(i+1);
    xx[i] = (hpt->GetBinContent(i+1,5)+hpt->GetBinContent(i+1,6))/2.;
    xxerr[i]=0;
  }
  TGraphErrors * g_20_30 = new TGraphErrors(numpt,xx,yy,xxerr,yyerr);
  g_20_30->SetMarkerStyle(marker_20_30);
  g_20_30->SetMarkerColor(color_20_30);
  g_20_30->SetLineColor(color_20_30);

  Int_t marker_30_40 = 25;
  Int_t color_30_40 = kViolet;
  TH1D * pos_30_40 = (TH1D *) pos->ProjectionY("pos_30_40",7,8);
  pos_30_40->Scale(1./resPos->GetBinContent(5));
  TH1D * neg_30_40 = (TH1D *) neg->ProjectionY("neg_30_40",7,8);
  neg_30_40->Scale(1./resNeg->GetBinContent(5));
  TH1D * v2_30_40 = (TH1D *) pos_30_40->Clone("v2_30_40");
  v2_30_40->Add(neg_30_40);
  v2_30_40->Scale(0.25);
  v2_30_40->SetMinimum(-0.1);
  v2_30_40->SetMaximum(0.4);
  v2_30_40->SetMarkerStyle(marker_30_40);
  v2_30_40->SetMarkerColor(color_30_40);
  v2_30_40->SetLineColor(color_30_40);
  for(int i = 0; i< numpt; i++ ) {
    yy[i] = v2_30_40->GetBinContent(i+1);
    yyerr[i] = v2_30_40->GetBinError(i+1);
    xx[i] = (hpt->GetBinContent(i+1,7)+hpt->GetBinContent(i+1,8))/2.;
    xxerr[i]=0;
  }
  TGraphErrors * g_30_40 = new TGraphErrors(numpt,xx,yy,xxerr,yyerr);
  g_30_40->SetMarkerStyle(marker_30_40);
  g_30_40->SetMarkerColor(color_30_40);
  g_30_40->SetLineColor(color_30_40);

  Int_t marker_40_50 = 20;
  Int_t color_40_50 = kSpring;
  TH1D * pos_40_50 = (TH1D *) pos->ProjectionY("pos_40_50",9,10);
  pos_40_50->Scale(1./resPos->GetBinContent(6));
  TH1D * neg_40_50 = (TH1D *) neg->ProjectionY("neg_40_50",9,10);
  neg_40_50->Scale(1./resNeg->GetBinContent(6));
  TH1D * v2_40_50 = (TH1D *) pos_40_50->Clone("v2_40_50");
  v2_40_50->Add(neg_40_50);
  v2_40_50->Scale(0.25);
  v2_40_50->SetMinimum(-0.1);
  v2_40_50->SetMaximum(0.4);
  v2_40_50->SetMarkerStyle(marker_40_50);
  v2_40_50->SetMarkerColor(color_40_50);
  v2_40_50->SetLineColor(color_40_50);
  for(int i = 0; i< numpt; i++ ) {
    yy[i] = v2_40_50->GetBinContent(i+1);
    yyerr[i] = v2_40_50->GetBinError(i+1);
    xx[i] = (hpt->GetBinContent(i+1,9)+hpt->GetBinContent(i+1,10))/2.;
    xxerr[i]=0;
  }
  TGraphErrors * g_40_50 = new TGraphErrors(numpt,xx,yy,xxerr,yyerr);
  g_40_50->SetMarkerStyle(marker_40_50);
  g_40_50->SetMarkerColor(color_40_50);
  g_40_50->SetLineColor(color_40_50);

  Int_t marker_50_60 = 29;
  Int_t color_50_60 = kOrange;
  TH1D * pos_50_60 = (TH1D *) pos->ProjectionY("pos_50_60",11,12);
  pos_50_60->Scale(1./resPos->GetBinContent(7));
  TH1D * neg_50_60 = (TH1D *) neg->ProjectionY("neg_50_60",11,12);
  neg_50_60->Scale(1./resNeg->GetBinContent(7));
  TH1D * v2_50_60 = (TH1D *) pos_50_60->Clone("v2_50_60");
  v2_50_60->Add(neg_50_60);
  v2_50_60->Scale(0.25);
  v2_50_60->SetMinimum(-0.1);
  v2_50_60->SetMaximum(0.4);
  v2_50_60->SetMarkerStyle(marker_50_60);
  v2_50_60->SetMarkerColor(color_50_60);
  v2_50_60->SetLineColor(color_50_60);
  for(int i = 0; i< numpt; i++ ) {
    yy[i] = v2_50_60->GetBinContent(i+1);
    yyerr[i] = v2_50_60->GetBinError(i+1);
    xx[i] = (hpt->GetBinContent(i+1,11)+hpt->GetBinContent(i+1,12))/2.;
    xxerr[i]=0;
    cout<<i<<" "<<xx[i]<<" "<<yy[i]<<endl;
  }
  TGraphErrors * g_50_60 = new TGraphErrors(numpt,xx,yy,xxerr,yyerr);
  g_50_60->SetMarkerStyle(marker_50_60);
  g_50_60->SetMarkerColor(color_50_60);
  g_50_60->SetLineColor(color_50_60);

  resPos->SetStats(kFALSE);
  resPos->SetTitle("Resolution Corrections");
  resPos->SetXTitle("Centrality");
  resPos->SetYTitle("ResCor");
  resPos->Draw();
  resNeg->SetLineColor(2);
  resNeg->Draw("same");
  TLegend * rleg = new TLegend(0.2,0.2,0.5,0.5,"EventPlane (3 sub-event)");
  rleg->AddEntry(resPos,"EvtPTracksPosEtaGap","lp");
  rleg->AddEntry(resNeg,"EvtPTracksNegEtaGap","lp");
  rleg->Draw();
  can->Print(Form("~/public_html/ResCor_%s.png",tag.Data()),"png");
  TCanvas * c2 = new TCanvas("c2","c2",800,600);
  gPad->SetGrid(1,1);
  TH1D * hfig = new TH1D("hfig",Form("|#eta|<1 (%s)",tag.Data()),100,0,6);
  hfig->SetMinimum(0.0);
  hfig->SetMaximum(0.25);
  hfig->SetStats(kFALSE);
  hfig->SetXTitle("p_{T} (GeV/c)");
  hfig->SetYTitle("v_{2}");
  hfig->Draw();
  g_0_5->Draw("p");
  g_5_10->Draw("p");
  g_10_20->Draw("p");
  g_20_30->Draw("p");
  g_30_40->Draw("p");
  g_40_50->Draw("p");
  g_50_60->Draw("p");

  TLegend * leg = new TLegend(0.65,0.65,0.85,0.88,"Centrality    (N_{part})");
  leg->AddEntry(v2_0_5,  Form("0-5     (%5.1f)",hNpartBin->GetBinContent(1)),"lp");
  leg->AddEntry(v2_5_10, Form("5-10    (%5.1f)",hNpartBin->GetBinContent(2)),"lp");
  leg->AddEntry(v2_10_20,Form("10-20   (%5.1f)",hNpartBin->GetBinContent(3)),"lp");
  leg->AddEntry(v2_20_30,Form("20-30   (%5.1f)",hNpartBin->GetBinContent(4)),"lp");
  leg->AddEntry(v2_30_40,Form("30-40   (%5.1f)",hNpartBin->GetBinContent(5)),"lp");
  leg->AddEntry(v2_40_50,Form("40-50   (%5.1f)",hNpartBin->GetBinContent(6)),"lp");
  leg->AddEntry(v2_50_60,Form("50-60   (%5.1f)",hNpartBin->GetBinContent(7)),"lp");
  leg->Draw();
  c2->Print(Form("~/public_html/v2_%s.png",tag.Data()),"png");

  cout<<"create c3"<<endl;
  TCanvas * c3 = new TCanvas("ptDist","ptDist",800,600);
  cout<<"canvas created"<<endl;
  TH1D * pt_0_5 = (TH1D *) dNdPt->ProjectionX("hpt_0_5",1,1);
  TH1D * hptframe = new TH1D("hptframe",tag.Data(),100,0,6);
  hptframe->SetMaximum(10000000);
  hptframe->SetMinimum(0.01);
  hptframe->SetStats(kFALSE);
  gPad->SetLogy();
  cout<<"SetLogy"<<endl;
  hptframe->SetXTitle("p_{T} (GeV/c)");
  hptframe->SetYTitle("#frac{1}{2#pi p_{t}}#frac{d^{2}N}{dp_{T}d#eta }");
  cout<<"Draw frame"<<endl;
  hptframe->Draw();
  cout<<"frame drawn"<<endl;
  int nptbins = pt_0_5->GetNbinsX();
  for(int i = 0; i<nptbins; i++ ) {
    if(xx[i]>0) {
      Double_t scale = 1./(2.*3.1415*xx[i]);
      yy[i]= scale*pt_0_5->GetBinContent(i+1);
      yyerr[i]=scale*pt_0_5->GetBinError(i+1);
    } else {
      yy[i] = 0;
      yyerr[i] = 0;
    }

  }
  TGraphErrors * gpt_0_5 = new TGraphErrors(nptbins,xx,yy,xxerr,yyerr);
  gpt_0_5->SetMarkerStyle(marker_0_5);
  gpt_0_5->SetMarkerColor(color_0_5);
  gpt_0_5->Draw("p");


  TH1D * pt_5_10 = (TH1D *) dNdPt->ProjectionX("hpt_5_10",2,2);
  for(int i = 0; i<nptbins; i++ ) {
    if(xx[i]>0) {
      Double_t scale = 1./(2.*3.1415*xx[i]);
      yy[i]= scale*pt_5_10->GetBinContent(i+1)/2;
      yyerr[i]=scale*pt_5_10->GetBinError(i+1)/2;
    } else { 
      yy[i] = 0;
      yyerr[i] = 0;
    }
    
  }
  TGraphErrors * gpt_5_10 = new TGraphErrors(nptbins,xx,yy,xxerr,yyerr);
  gpt_5_10->SetMarkerStyle(marker_5_10);
  gpt_5_10->SetMarkerColor(color_5_10);
  gpt_5_10->Draw("p");

  TH1D * pt_10_20 = (TH1D *) dNdPt->ProjectionX("hpt_10_20",3,4);
  pt_10_20->Scale(0.5);
  for(int i = 0; i<nptbins; i++ ) {
    if(xx[i] > 0) {
      Double_t scale = 1./(2.*3.1415*xx[i]);
      yy[i]= scale*pt_10_20->GetBinContent(i+1)/2/2;
      yyerr[i]=scale*pt_10_20->GetBinError(i+1)/2/2;
    } else { 
      yy[i] = 0;
      yyerr[i] = 0;
    }
  }
  TGraphErrors * gpt_10_20 = new TGraphErrors(nptbins,xx,yy,xxerr,yyerr);
  gpt_10_20->SetMarkerStyle(marker_10_20);
  gpt_10_20->SetMarkerColor(color_10_20);
  gpt_10_20->Draw("p");

  TH1D * pt_20_30 = (TH1D *) dNdPt->ProjectionX("hpt_20_30",5,6);
  pt_20_30->Scale(0.5);
  for(int i = 0; i<nptbins; i++ ) {
    if( xx[i] > 0 ) {
      Double_t scale = 1./(2.*3.1415*xx[i]);
      yy[i]= scale*pt_20_30->GetBinContent(i+1)/2/2/2;
      yyerr[i]=scale*pt_20_30->GetBinError(i+1)/2/2/2;
    } else {
      yy[i] = 0;
      yyerr[i] = 0;
    }
  }
  TGraphErrors * gpt_20_30 = new TGraphErrors(nptbins,xx,yy,xxerr,yyerr);
  gpt_20_30->SetMarkerStyle(marker_20_30);
  gpt_20_30->SetMarkerColor(color_20_30);
  gpt_20_30->Draw("p");

  TH1D * pt_30_40 = (TH1D *) dNdPt->ProjectionX("hpt_30_40",7,8);
  pt_30_40->Scale(0.5);
  for(int i = 0; i<nptbins; i++ ) {
    if(xx[i]>0) {
      Double_t scale = 1./(2.*3.1415*xx[i]);
      yy[i]= scale*pt_30_40->GetBinContent(i+1)/2/2/2/2;
      yyerr[i]=scale*pt_30_40->GetBinError(i+1)/2/2/2/2;
    } else {
      yy[i] = 0;
      yyerr[i] = 0;
    }
  }
  TGraphErrors * gpt_30_40 = new TGraphErrors(nptbins,xx,yy,xxerr,yyerr);
  gpt_30_40->SetMarkerStyle(marker_30_40);
  gpt_30_40->SetMarkerColor(color_30_40);
  gpt_30_40->Draw("p");

  TH1D * pt_40_50 = (TH1D *) dNdPt->ProjectionX("hpt_40_50",9,10);
  pt_40_50->Scale(0.5);
  for(int i = 0; i<nptbins; i++ ) {
    if( xx[i] ) {
      Double_t scale = 1./(2.*3.1415*xx[i]);
      yy[i]= scale*pt_40_50->GetBinContent(i+1)/2/2/2/2/2;
      yyerr[i]=scale*pt_40_50->GetBinError(i+1)/2/2/2/2/2;
    } else { 
      yy[i] = 0;
      yyerr[i] = 0;
    }
  }
  TGraphErrors * gpt_40_50 = new TGraphErrors(nptbins,xx,yy,xxerr,yyerr);
  gpt_40_50->SetMarkerStyle(marker_40_50);
  gpt_40_50->SetMarkerColor(color_40_50);
  gpt_40_50->Draw("p");

  TH1D * pt_50_60 = (TH1D *) dNdPt->ProjectionX("hpt_50_60",11,12);
  pt_50_60->Scale(0.5);
  for(int i = 0; i<nptbins; i++ ) {
    if(xx[i]>0) {
      Double_t scale = 1./(2.*3.1415*xx[i]);
      yy[i]= scale*pt_50_60->GetBinContent(i+1)/2/2/2/2/2/2;
      yyerr[i]=scale*pt_50_60->GetBinError(i+1)/2/2/2/2/2/2;
    } else {
      yy[i] = 0;
      yyerr[i] = 0;
   }
  }
  TGraphErrors * gpt_50_60 = new TGraphErrors(nptbins,xx,yy,xxerr,yyerr);
  gpt_50_60->SetMarkerStyle(marker_50_60);
  gpt_50_60->SetMarkerColor(color_50_60);
  gpt_50_60->Draw("p");


  c3->Print(Form("~/public_html/ptDist_%s.png",tag.Data()),"png");

  TLegend * leg2 = new TLegend(0.60,0.50,0.85,0.88,"Centrality    (N_{part})");
  leg2->AddEntry(v2_0_5,  Form("0-5     (%5.1f)",hNpartBin->GetBinContent(1)),"lp");
  leg2->AddEntry(v2_5_10, Form("5-10    (%5.1f)",hNpartBin->GetBinContent(2)),"lp");
  leg2->AddEntry(v2_10_20,Form("10-20   (%5.1f)",hNpartBin->GetBinContent(3)),"lp");
  leg2->AddEntry(v2_20_30,Form("20-30   (%5.1f)",hNpartBin->GetBinContent(4)),"lp");
  leg2->AddEntry(v2_30_40,Form("30-40   (%5.1f)",hNpartBin->GetBinContent(5)),"lp");
  leg2->AddEntry(v2_40_50,Form("40-50   (%5.1f)",hNpartBin->GetBinContent(6)),"lp");
  leg2->AddEntry(v2_50_60,Form("50-60   (%5.1f)",hNpartBin->GetBinContent(7)),"lp");
  leg2->Draw();

  std::ofstream file;
  file.open("v2Results.txt");
  file<<tag.Data()<<endl<<endl;
  file<<Form("0-5     (Npart = %5.1f)",hNpartBin->GetBinContent(1))<<endl;
  for(int i = 0; i< nptbins;  i++ ) {
    file<<setprecision(3)<<xx[i]<<"\t"<<setprecision(3)<<v2_0_5->GetBinContent(i+1)<<"\t"<<v2_0_5->GetBinError(i+1)<<endl;
  }
  file<<endl;

  file<<Form("5-10     (Npart = %5.1f)",hNpartBin->GetBinContent(2))<<endl;
  for(int i = 0; i< nptbins;  i++ ) {
    file<<setprecision(3)<<xx[i]<<"\t"<<setprecision(3)<<v2_5_10->GetBinContent(i+1)<<"\t"<<v2_5_10->GetBinError(i+1)<<endl;
  }
  file<<endl;

  file<<Form("10-20     (Npart = %5.1f)",hNpartBin->GetBinContent(3))<<endl;
  for(int i = 0; i< nptbins;  i++ ) {
    file<<setprecision(3)<<xx[i]<<"\t"<<setprecision(3)<<v2_10_20->GetBinContent(i+1)<<"\t"<<v2_10_20->GetBinError(i+1)<<endl;
  }
  file<<endl;

  file<<Form("20-30     (Npart = %5.1f)",hNpartBin->GetBinContent(4))<<endl;
  for(int i = 0; i< nptbins;  i++ ) {
    file<<setprecision(3)<<xx[i]<<"\t"<<setprecision(3)<<v2_20_30->GetBinContent(i+1)<<"\t"<<v2_20_30->GetBinError(i+1)<<endl;
  }
  file<<endl;

  file<<Form("30-40     (Npart = %5.1f)",hNpartBin->GetBinContent(5))<<endl;
  for(int i = 0; i< nptbins;  i++ ) {
    file<<setprecision(3)<<xx[i]<<"\t"<<setprecision(3)<<v2_30_40->GetBinContent(i+1)<<"\t"<<v2_30_40->GetBinError(i+1)<<endl;
  }
  file<<endl;

  file<<Form("40-50     (Npart = %5.1f)",hNpartBin->GetBinContent(6))<<endl;
  for(int i = 0; i< nptbins;  i++ ) {
    file<<setprecision(3)<<xx[i]<<"\t"<<setprecision(3)<<v2_40_50->GetBinContent(i+1)<<"\t"<<v2_40_50->GetBinError(i+1)<<endl;
  }
  file<<endl;

  file<<Form("50-60     (Npart = %5.1f)",hNpartBin->GetBinContent(7))<<endl;
  for(int i = 0; i< nptbins;  i++ ) {
    file<<setprecision(3)<<xx[i]<<"\t"<<setprecision(3)<<v2_50_60->GetBinContent(i+1)<<"\t"<<v2_50_60->GetBinError(i+1)<<endl;
  }
  file<<endl;
  file.close();
}
