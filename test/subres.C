TH2D * EPCorr[20];
TGraph * MidEta;
TGraph * MidEtaGen;
TGraph * NegEta;
TGraph * NegEtaGen;
TGraph * PosEta;
TGraph * PosEtaGen;
TGraph * v2EtaGraph[8];
#include "/net/hisrv0001/home/sanders/CMSSW_3_7_0/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
static const Double_t EtaGap = 0.4; //full gap between subevents
static const Double_t AutoEtaGap = 0.1;  // full gap of autocorrelation rejected particles
static const Double_t TrackEtaWidth = 8.;
static const Double_t HFEtaWidth = 10.;
double multcorrTrack = (TrackEtaWidth-AutoEtaGap)/TrackEtaWidth;
double multcorrHF = (TrackEtaWidth - AutoEtaGap)/TrackEtaWidth;
void subres(){
  TString tag = "SubRes_eta4_4cuts_gap1_mult";
  TFile * tf = new TFile("/net/hibat0007/d00/scratch/sanders/370/rpflat_combined.root");
  TF1 * f2 = new TF1("f2","pol5",0,3);
  f2->SetParameters(0.,0.626657,0.,-0.09694,0.02754,-0.002283);
  TF1 * if2 = new TF1("if2","pol7",0,1);
  if2->SetParameters(-1.588838e-02,2.783860e+00,-1.966387e+01,1.314961e+02,-4.271861e+02,7.267611e+02,-6.191080e+02,2.094568e+02);
  TList * list = ((TDirectory *) tf->Get("hiEvtPlaneFlat"))->GetListOfKeys();
  int indx = 0;
  TH1D * subr[50];
  TH1D * genr[50];
  Double_t subrval[50];
  Double_t subrerr[50];
  Double_t genrval[50];
  Double_t genrerr[50];
  Double_t fullrval[50];
  Double_t fullrerr[50];
  TH1D * fullmult;
  TH1D * sub1mult;
  TH1D * sub2mult;

  TGraphErrors * res[50];
  TGraphErrors * resR[50];
  Double_t xval[50];
  Double_t yval[50];
  Double_t xerr[50];
  Double_t yerr[50];
  TString names[50];
  Int_t indx = 0;
  TCanvas * can = new TCanvas("SubEvtRes","SubEvtRes",1300,800);
  can->Divide(3,2);
  can->cd(6);
  TLatex * t = new TLatex(.2,.6,Form("Subevent #eta gap of #pm %4.2f",EtaGap/2.));
  t->Draw();
  TLatex * t2 = new TLatex(.2,.4,Form("Autocorrelation #eta gap of #pm %4.2f",AutoEtaGap/2.));
  t2->Draw();
  can->cd(1);
  gPad->SetGrid(1,1);
  TH1D * h = new TH1D("h","ResCor(random)/True ResCor",500,0,100);
  h->SetStats(kFALSE);
  h->SetXTitle("Centrality");
  h->SetYTitle("Resolution Correction");
  h->SetMaximum(1.6);
  h->SetMinimum(0.4);
  h->Draw();

  Int_t idx=0;
  while(indx>=0) {
    TString name = list->At(indx)->GetName();
    if(name.Contains("1")) goto next;
    if(!name.Contains("Track") && !name.Contains("etCaloHF")) goto next;
    if(name.Contains("Ch")) goto next;
    if(name.Contains("CaloHFP")) goto next;
    if(name.Contains("CaloHFM")) goto next;
    if(name.Contains("EvtPlaneFromTracksPosEta")) goto next;
    if(name.Contains("EvtPlaneFromTracksNegEta")) goto next;
    
    //if(name.Contains("etC")) goto next;
    
    if(name.Contains("cent")) break;
    if(!name.Contains("MidEtaTrackRescor")) {
      cout<<name.Data()<<endl;
      names[idx] = name;
      Int_t npts = 0;
      for(int i = 0; i<50; i++) {
	xval[i]=0; 
	yval[i]=0;
	yerr[i]=0;
	subrval[i]=0;
	subrerr[i]=0;
	genrval[i]=0;
	genrerr[i]=0;
	fullrval[i]=0;
	fullrerr[i]=0;
      }
      sub1mult = (TH1D *) tf->Get(Form("v2analyzer/EventPlanes/%s/Mult/Mult1",name.Data()));
      sub2mult = (TH1D *) tf->Get(Form("v2analyzer/EventPlanes/%s/Mult/Mult2",name.Data()));
      fullmult = (TH1D *) tf->Get(Form("v2analyzer/EventPlanes/%s/Mult/Mult",name.Data()));
      sub1mult->Divide((TH1D *) tf->Get(Form("v2analyzer/EventPlanes/%s/Mult/Mult1Cnt",name.Data())));
      sub2mult->Divide((TH1D *) tf->Get(Form("v2analyzer/EventPlanes/%s/Mult/Mult2Cnt",name.Data())));
      fullmult->Divide((TH1D *) tf->Get(Form("v2analyzer/EventPlanes/%s/Mult/MultCnt",name.Data())));
      for (int i = 0; i<15; i++) {
	subr[i] = (TH1D *) tf->Get(Form("v2analyzer/EventPlanes/%s/SubRes/SubRes_%d",name.Data(),i));
	genr[i] = (TH1D *) tf->Get(Form("v2analyzer/EventPlanes/%s/GenRes/GenRes_%d",name.Data(),i));
	double multcorr = 1.;
	if(name.Contains("EvtPlaneFromTracksEta")) multcorr = (TrackEtaWidth-AutoEtaGap)/TrackEtaWidth;
	if(name.Contains("etCaloHF"))  multcorr = (TrackEtaWidth - AutoEtaGap)/TrackEtaWidth;
	Double_t ratio = multcorr*fullmult->GetBinContent(i)/(2.*(sub1mult->GetBinContent(i)*sub2mult->GetBinContent(i))/(sub1mult->GetBinContent(i)+sub2mult->GetBinContent(i)));
	subrval[i]=subr[i]->GetMean();
	subrerr[i]=subr[i]->GetMeanError();
	genrval[i]=genr[i]->GetMean();
	genrerr[i]=genr[i]->GetMeanError();
	if(subrval[i]>0 && subrerr[i]>0) {
	  Double_t chi = if2->Eval(sqrt(subrval[i]));
	  fullrval[i] = f2->Eval(sqrt(ratio)*chi);
	  fullrerr[i] = 0.5*fullrval[i]*subrerr[i]/subrval[i];
	  xval[npts] = 5*i+2.5;
	  xerr[npts] = 2.5;
	  yval[npts] = fullrval[i]/genrval[i];
	  yerr[npts] = fullrerr[i]/genrval[i];
	  ++npts;
	}
      }
      res[idx] = new TGraphErrors(npts,xval,yval,xerr,yerr);
      res[idx]->SetName(name.Data());
      res[idx]->SetTitle(name.Data());
      res[idx]->SetMarkerStyle(29);
      resR[idx] = new TGraphErrors(npts,xval,fullrval,xerr,fullrerr);
      resR[idx]->SetName(name.Data());
      resR[idx]->SetTitle(name.Data());
      resR[idx]->SetMarkerStyle(29);
      if(name.Contains("Track") ) {
	if(name.Contains("EvtPlaneFromTracksMidEta")) {
	  res[idx]->SetMarkerStyle(20); 
	  res[idx]->SetMarkerSize(1.2);
	  res[idx]->SetMarkerColor(kBlue);
	  res[idx]->SetLineColor(kBlue);
	  resR[idx]->SetMarkerStyle(20); 
	  resR[idx]->SetMarkerSize(1.2);
	  resR[idx]->SetMarkerColor(kBlue);
	  resR[idx]->SetLineColor(kBlue);
	}
	if(name.Contains("EvtPTracksPosEtaGap")) {
	  res[idx]->SetMarkerStyle(21); 
	  res[idx]->SetMarkerSize(1.2);
	  res[idx]->SetMarkerColor(kRed);
	  res[idx]->SetLineColor(kRed);
	  resR[idx]->SetMarkerStyle(21); 
	  resR[idx]->SetMarkerSize(1.2);
	  resR[idx]->SetMarkerColor(kRed);
	  resR[idx]->SetLineColor(kRed);
	}
	if(name.Contains("EvtPTracksNegEtaGap")) {
	  res[idx]->SetMarkerStyle(22); 
	  res[idx]->SetMarkerSize(1.2);
	  res[idx]->SetMarkerColor(kGreen);
	  res[idx]->SetLineColor(kGreen);
	  resR[idx]->SetMarkerStyle(22); 
	  resR[idx]->SetMarkerSize(1.2);
	  resR[idx]->SetMarkerColor(kGreen);
	  resR[idx]->SetLineColor(kGreen);
	}
	if(name.Contains("EvtPlaneFromTracksEta")) {
	  res[idx]->SetMarkerStyle(26); 
	  res[idx]->SetMarkerSize(1.2);
	  res[idx]->SetMarkerColor(kMagenta);
	  res[idx]->SetLineColor(kMagenta);
	  resR[idx]->SetMarkerStyle(26); 
	  resR[idx]->SetMarkerSize(1.2);
	  resR[idx]->SetMarkerColor(kMagenta);
	  resR[idx]->SetLineColor(kMagenta);
	}
      }
      if(name.Contains("Calo") ) {
	  res[idx]->SetMarkerStyle(27); 
	  res[idx]->SetMarkerSize(1.2);
	  res[idx]->SetMarkerColor(kRed);
	  res[idx]->SetLineColor(kRed);
	  resR[idx]->SetMarkerStyle(27); 
	  resR[idx]->SetMarkerSize(1.2);
	  resR[idx]->SetMarkerColor(kRed);
	  resR[idx]->SetLineColor(kRed);
      }
      if(name.Contains("Ecal") ) {
	res[idx]->SetMarkerColor(kCyan);
	res[idx]->SetLineColor(kCyan);
      } 
      if(name.Contains("Hcal") ) {
	res[idx]->SetMarkerColor(kGreen);
	res[idx]->SetLineColor(kGreen);
      } 
      if(name.Contains("1") ) {
	res[idx]->SetMarkerColor(kMagenta);
	res[idx]->SetLineColor(kMagenta);
      } 
     res[idx]->Draw("P same");
     ++idx;
    }
  next:
    if(list->At(indx)==list->Last())
      indx = -1;
    else
      ++indx;
  }
  can->cd(5);
  TLegend * leg1 = new TLegend(0.1,0.20,0.5,0.8);
  for (int i = 0; i< idx; i++) {
    leg1->AddEntry(res[i],names[i].Data(),"p");
  } 
  leg1->Draw();
  
  can->cd(2);
  gPad->SetGrid(1,1);
  TH1D * v2Gen[8];
  TH2D * v2GenCos = tf->Get("v2analyzer/v2/v2Reco/v2GenCos");
  TH2D * v2GenCnt = tf->Get("v2analyzer/v2/v2Reco/v2GenCnt");
  v2GenCos->Divide(v2GenCnt);
  
  
  TH1D * etCaloHFCos[8];
  TH1D * etCaloHFCnt[8];
  for( int i = 0; i<8; i++) {
    v2Gen[i] = v2GenCos->ProjectionX(Form("v2Gen_%d",i),i+1,i+1);
    TH2D * hcos = (TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/etCaloHF/cos_etCaloHF_%d",i));
    TH2D * hcnt = (TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/etCaloHF/cnt_etCaloHF_%d",i));
    etCaloHFCos[i] = (TH1D *) (hcos->ProjectionX(Form("etCaloHFCos_%d",i))); 
    etCaloHFCnt[i] = (TH1D *) (hcnt->ProjectionX(Form("etCaloHFCnt_%d",i))); 
    etCaloHFCos[i]->Divide(etCaloHFCnt[i]);
    etCaloHFCos[i]->SetLineColor(i+1);
    etCaloHFCos[i]->Divide(v2Gen[i]);
    for(int j = 0; j<20; j++) {
      if(fullrval[j]>0) etCaloHFCos[i]->SetBinContent(j+1,etCaloHFCos[i]->GetBinContent(j+1)/fullrval[j]);
    }
    if(i==0) {
      etCaloHFCos[i]->SetMinimum(0.5);
      etCaloHFCos[i]->SetMaximum(1.5);
      etCaloHFCos[i]->SetTitle("v2_{etCaloHF}/v2_{Gen}");
      etCaloHFCos[i]->Draw();
      
    } else {
      etCaloHFCos[i]->Draw("same");
    } 
  }
  
  can->cd(3);
  gPad->SetGrid(1,1);
  TH1D * TrackCos[8];
  TH1D * TrackCnt[8];
  for( int i = 0; i<8; i++) {
    TH2D * hcos = (TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/EvtPlaneFromTracksEta/cos_EvtPlaneFromTracksEta_%d",i));
    TH2D * hcnt = (TH2D *) tf->Get(Form("v2analyzer/v2/v2Reco/EvtPlaneFromTracksEta/cnt_EvtPlaneFromTracksEta_%d",i));
    TrackCos[i] = (TH1D *) (hcos->ProjectionX(Form("EvtPlaneFromTracksEtaCos_%d",i))); 
    TrackCnt[i] = (TH1D *) (hcnt->ProjectionX(Form("EvtPlaneFromTracksEtaCnt_%d",i))); 
    TrackCos[i]->Divide(etCaloHFCnt[i]);
    TrackCos[i]->SetLineColor(i+1);
    TrackCos[i]->Divide(v2Gen[i]);
    for(int j = 0; j<20; j++ ) {
      if(fullrval[j]>0) TrackCos[i]->SetBinContent(j+1,TrackCos[i]->GetBinContent(j+1)/fullrval[j]);
    }
    if(i==0) {
      TrackCos[0]->SetMinimum(0.5);
      TrackCos[0]->SetMaximum(1.5);
      TrackCos[i]->SetTitle("v2_{EvtPlaneFromTracksEta}/v2_{Gen}");
      TrackCos[i]->Draw();
    } else {
      TrackCos[i]->Draw("same");
    }
  }
  can->cd(5);
  TLegend * leg2 = new TLegend(0.55,0.20,0.9,0.8);
  for (int i = 0; i< 8; i++) {
    leg2->AddEntry(TrackCos[i],Form("%4.1f #leq #eta #leq %4.1f",-2+0.5*i,-1.5+0.5*i),"pl");
  } 
  leg2->Draw();
  
  can->cd(4);
  gPad->SetGrid(1,1);
  h3 = (TH1D *) h->Clone("h3");
  h3->SetMinimum(0.0);
  h3->SetMaximum(1.0);
  h3->SetTitle("ResCor from subevent");
  h3->Draw();
  for (int i = 0; i< idx; i++) {
    resR[i]->Draw("p same");
  } 
  can->Print(Form("%s.png",tag.Data()),"png");
}
