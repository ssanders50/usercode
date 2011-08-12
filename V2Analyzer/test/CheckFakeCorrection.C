#include <iomanip>
#include <fstream>
TString TagLabel;
TGraphErrors * InputType1(TString fileName, Int_t & Num);
TGraphErrors * InputType1(int  minCent, int maxCent, double eta, TString tag, TString EP1, TString EP2,Int_t & Num) {
  TString name = Form("results/v2pt_EP_%d-%d_%02d-%s_%s_%s.txt",minCent,maxCent,(int)(10*eta),tag.Data(),EP1.Data(),EP2.Data());
  if(EP2.Length()==0) name = Form("results/v2pt_EP_%d-%d_%02d_%s_%s.txt",minCent,maxCent,(int)(10*eta),tag.Data(),EP1.Data());
  TagLabel = Form("%s_%s_%s",tag.Data(),EP1.Data(),EP2.Data());
  return InputType1(name,Num);
}
void CheckFakeCorrection(){

  int minc = 5;
  int maxc =10;

  Double_t lmargin1 = 0.1;
  Double_t rmargin1 = 0.01;
  Double_t lmargin2 = 0.1;
  Double_t rmargin2 = 0.05;
  Double_t xaxislen = (1.-(rmargin1+lmargin1+rmargin2+lmargin2))/2.;
  Double_t xpad1 = lmargin1+xaxislen+rmargin1;
  Double_t xpad2 = 1-xpad1;
  TCanvas * can = new TCanvas("FakeCorrection","FakeCorrection",1300,750);
  can->cd();
  TPad * pad1 = new TPad("pad1","pad1",0.0,0.0,xpad1,1.0);
  pad1->Draw();
  pad1->SetLeftMargin(lmargin1/xpad1);
  pad1->SetRightMargin(rmargin1/xpad1);
  pad1->cd();
  TH1D * h = new TH1D("h","",500,0,1.2);
  h->SetMaximum(0.07);
  h->SetMinimum(0.0);
  h->SetXTitle("p_{T}");
  h->SetYTitle("v_{2}{EP}");
  h->SetStats(kFALSE);
  h->Draw();
  //TString tag = "Asym-gt10dz40chi-minpt2";
  TString tag = "Asym-gt10dz40chi12mean";
  TString ep1 = "etHF";
  TString ep2 = "";
  TGraphErrors * dz10_0_5    = InputType1(0,    5,0.8, tag,  ep1, ep2,    13);
  TGraphErrors * dz10_5_10   = InputType1(5,   10,0.8, tag,  ep1, ep2,    13);
  TGraphErrors * dz10_10_15  = InputType1(10,  15,0.8, tag,  ep1, ep2,    13);
  //tag = "Asym-gt6dz20chi-minpt2";
  tag = "Asym-gt6dz20chi12mean";
  TGraphErrors * dz3_0_5 = InputType1(0, 5, 0.8, tag, ep1, ep2, 13);
  TGraphErrors * dz3_5_10 = InputType1(5, 10, 0.8, tag, ep1, ep2, 13);
  TGraphErrors * dz3_10_15  = InputType1(10,  15,0.8, tag,  ep1, ep2,    13);

  dz10_0_5->SetMarkerStyle(20);
  dz10_0_5->SetMarkerColor(kBlue);
  dz3_0_5->SetMarkerStyle(24);
  dz3_0_5->SetMarkerColor(kBlue);
  dz10_5_10->SetMarkerStyle(21);
  dz10_5_10->SetMarkerColor(kRed);
  dz3_5_10->SetMarkerStyle(25);
  dz3_5_10->SetMarkerColor(kRed);
  dz10_10_15->SetMarkerStyle(29);
  dz10_10_15->SetMarkerColor(kMagenta+1);
  dz3_10_15->SetMarkerStyle(30);
  dz3_10_15->SetMarkerColor(kMagenta+1);
  dz10_0_5->Draw("p");
  dz3_0_5->Draw("p");
  dz10_5_10->Draw("p");
  dz3_5_10->Draw("p");
  dz10_10_15->Draw("p");
  dz3_10_15->Draw("p");
  TLegend * leg = new TLegend(0.24,0.6,0.44,0.9);
  leg->AddEntry(dz10_0_5,  Form("0-5%c (10/40)",'%'),"p");
  leg->AddEntry(dz10_5_10, Form("5-10%c (10/40)",'%'),"p");
  leg->AddEntry(dz10_10_15,Form("10-15%c (10/40)",'%'),"p");
  leg->AddEntry(dz3_0_5,   Form("0-5%c (3/20)",'%'),"p");
  leg->AddEntry(dz3_5_10,  Form("5-10%c (3/20)",'%'),"p");
  leg->AddEntry(dz3_10_15, Form("10-15%c (3/20)",'%'),"p");
  leg->Draw();
  TLatex * text = new TLatex(0.7,0.012,"Uncorrected");
  text->Draw();
  gPad->SetGrid(1,1);
  //tag = "Asym-gt10dz40chi-minpt2";
  tag = "Asym-gt10dz40chi12mean";
  ep1 = "etHF_FAKECORRECTED";
  TGraphErrors * dz10f_0_5    = InputType1(0,  5,0.8, tag,  ep1, ep2,    13);
  TGraphErrors * dz10f_5_10   = InputType1(5,  10,0.8, tag,  ep1, ep2,    13);
  TGraphErrors * dz10f_10_15  = InputType1(10,  15,0.8, tag,  ep1, ep2,    13);
  //tag = "Asym-gt6dz20chi-minpt2";
  tag = "Asym-gt6dz20chi12mean";
  TGraphErrors * dz3f_0_5   = InputType1(0, 5, 0.8, tag, ep1, ep2, 13);
  TGraphErrors * dz3f_5_10  = InputType1(5, 10, 0.8, tag, ep1, ep2, 13);
  TGraphErrors * dz3f_10_15 = InputType1(10, 15, 0.8, tag, ep1, ep2, 13);

  Double_t * x_dz10f_0_5     = dz10_0_5->GetX();
  Double_t * y_dz10f_0_5     = dz10_0_5->GetY();
  Double_t * err_dz10f_0_5   = dz10_0_5->GetEY();
  Double_t * x_dz10f_5_10    = dz10_5_10->GetX();
  Double_t * y_dz10f_5_10    = dz10_5_10->GetY();
  Double_t * err_dz10f_5_10  = dz10_5_10->GetEY();
  Double_t * x_dz10f_10_15   = dz10_10_15->GetX();
  Double_t * y_dz10f_10_15   = dz10_10_15->GetY();
  Double_t * err_dz10f_10_15 = dz10_10_15->GetEY();

  Double_t * x_dz3f_0_5     = dz3_0_5->GetX();
  Double_t * y_dz3f_0_5     = dz3_0_5->GetY();
  Double_t * err_dz3f_0_5   = dz3_0_5->GetEY();
  Double_t * x_dz3f_5_10    = dz3_5_10->GetX();
  Double_t * y_dz3f_5_10    = dz3_5_10->GetY();
  Double_t * err_dz3f_5_10  = dz3_5_10->GetEY();
  Double_t * x_dz3f_10_15   = dz3_10_15->GetX();
  Double_t * y_dz3f_10_15   = dz3_10_15->GetY();
  Double_t * err_dz3f_10_15 = dz3_10_15->GetEY();



  dz10f_0_5->SetMarkerStyle(20);
  dz10f_0_5->SetMarkerColor(kBlue);
  dz3f_0_5->SetMarkerStyle(24);
  dz3f_0_5->SetMarkerColor(kBlue);

  dz10f_5_10->SetMarkerStyle(21);
  dz10f_5_10->SetMarkerColor(kRed);
  dz3f_5_10->SetMarkerStyle(25);
  dz3f_5_10->SetMarkerColor(kRed);

  dz10f_10_15->SetMarkerStyle(29);
  dz10f_10_15->SetMarkerColor(kMagenta+1);
  dz3f_10_15->SetMarkerStyle(30);
  dz3f_10_15->SetMarkerColor(kMagenta+1);
  can->cd();
  TPad * pad2 = new TPad("pad2","pad2",xpad2,0.,1.,1.0);
  pad2->Draw();
  pad2->SetLeftMargin(lmargin2/xpad2);
  pad2->SetRightMargin(rmargin2/xpad2);
  pad2->cd();
  h->Draw();
  dz10f_0_5->Draw("p");
  dz3f_0_5->Draw("p");
  dz10f_5_10->Draw("p");
  dz3f_5_10->Draw("p");
  dz10f_10_15->Draw("p");
  dz3f_10_15->Draw("p");

  pad2->SetGrid(1,1);
  TLatex * text2 = new TLatex(0.7,0.012,"Corrected");
  text2->Draw();
  can->Print("FakeCorrection.png","png");




  Double_t * x_dz10_0_5     = dz10_0_5->GetX();
  Double_t * y_dz10_0_5     = dz10_0_5->GetY();
  Double_t * err_dz10_0_5   = dz10_0_5->GetEY();
  Double_t * x_dz10_5_10    = dz10_5_10->GetX();
  Double_t * y_dz10_5_10    = dz10_5_10->GetY();
  Double_t * err_dz10_5_10  = dz10_5_10->GetEY();
  Double_t * x_dz10_10_15   = dz10_10_15->GetX();
  Double_t * y_dz10_10_15   = dz10_10_15->GetY();
  Double_t * err_dz10_10_15 = dz10_10_15->GetEY();

  Double_t * x_dz3_0_5     = dz3_0_5->GetX();
  Double_t * y_dz3_0_5     = dz3_0_5->GetY();
  Double_t * err_dz3_0_5   = dz3_0_5->GetEY();
  Double_t * x_dz3_5_10    = dz3_5_10->GetX();
  Double_t * y_dz3_5_10    = dz3_5_10->GetY();
  Double_t * err_dz3_5_10  = dz3_5_10->GetEY();
  Double_t * x_dz3_10_15   = dz3_10_15->GetX();
  Double_t * y_dz3_10_15   = dz3_10_15->GetY();
  Double_t * err_dz3_10_15 = dz3_10_15->GetEY();

  for(Int_t i = 0; i<dz10_0_5->GetN();i++) {
    err_dz3_0_5[i]=y_dz3_0_5[i]*TMath::Sqrt( pow(err_dz10_0_5[i]/y_dz10_0_5[i],2)+pow(err_dz3_0_5[i]/y_dz3_0_5[i],2));
    y_dz3_0_5[i]/=y_dz10_0_5[i];

    err_dz3_5_10[i]=y_dz3_5_10[i]*TMath::Sqrt( pow(err_dz10_5_10[i]/y_dz10_5_10[i],2)+pow(err_dz3_5_10[i]/y_dz3_5_10[i],2));
    y_dz3_5_10[i]/=y_dz10_5_10[i];

    err_dz3_10_15[i]=y_dz3_10_15[i]*TMath::Sqrt( pow(err_dz10_10_15[i]/y_dz10_10_15[i],2)+pow(err_dz3_10_15[i]/y_dz3_10_15[i],2));
    y_dz3_10_15[i]/=y_dz10_10_15[i];
  }
  TGraphErrors * r_dz3_0_5 = new TGraphErrors(dz10_0_5->GetN(),x_dz3_0_5,y_dz3_0_5,0,err_dz3_0_5);
  TGraphErrors * r_dz3_5_10 = new TGraphErrors(dz10_0_5->GetN(),x_dz3_5_10,y_dz3_5_10,0,err_dz3_5_10);
  TGraphErrors * r_dz3_10_15 = new TGraphErrors(dz10_0_5->GetN(),x_dz3_10_15,y_dz3_10_15,0,err_dz3_10_15);
  r_dz3_0_5->SetMarkerStyle(24);
  r_dz3_0_5->SetMarkerColor(kBlue);

  r_dz3_5_10->SetMarkerStyle(25);
  r_dz3_5_10->SetMarkerColor(kRed);

  r_dz3_10_15->SetMarkerStyle(30);
  r_dz3_10_15->SetMarkerColor(kMagenta+1);


  Double_t * x_dz10f_0_5     = dz10f_0_5->GetX();
  Double_t * y_dz10f_0_5     = dz10f_0_5->GetY();
  Double_t * err_dz10f_0_5   = dz10f_0_5->GetEY();
  Double_t * x_dz10f_5_10    = dz10f_5_10->GetX();
  Double_t * y_dz10f_5_10    = dz10f_5_10->GetY();
  Double_t * err_dz10f_5_10  = dz10f_5_10->GetEY();
  Double_t * x_dz10f_10_15   = dz10f_10_15->GetX();
  Double_t * y_dz10f_10_15   = dz10f_10_15->GetY();
  Double_t * err_dz10f_10_15 = dz10f_10_15->GetEY();

  Double_t * x_dz3f_0_5     = dz3f_0_5->GetX();
  Double_t * y_dz3f_0_5     = dz3f_0_5->GetY();
  Double_t * err_dz3f_0_5   = dz3f_0_5->GetEY();
  Double_t * x_dz3f_5_10    = dz3f_5_10->GetX();
  Double_t * y_dz3f_5_10    = dz3f_5_10->GetY();
  Double_t * err_dz3f_5_10  = dz3f_5_10->GetEY();
  Double_t * x_dz3f_10_15   = dz3f_10_15->GetX();
  Double_t * y_dz3f_10_15   = dz3f_10_15->GetY();
  Double_t * err_dz3f_10_15 = dz3f_10_15->GetEY();

  for(Int_t i = 0; i<dz10f_0_5->GetN();i++) {
    err_dz3f_0_5[i]=y_dz3f_0_5[i]*TMath::Sqrt( pow(err_dz10f_0_5[i]/y_dz10f_0_5[i],2)+pow(err_dz3f_0_5[i]/y_dz3f_0_5[i],2));
    y_dz3f_0_5[i]/=y_dz10f_0_5[i];

    err_dz3f_5_10[i]=y_dz3f_5_10[i]*TMath::Sqrt( pow(err_dz10f_5_10[i]/y_dz10f_5_10[i],2)+pow(err_dz3f_5_10[i]/y_dz3f_5_10[i],2));
    y_dz3f_5_10[i]/=y_dz10f_5_10[i];

    err_dz3f_10_15[i]=y_dz3f_10_15[i]*TMath::Sqrt( pow(err_dz10f_10_15[i]/y_dz10f_10_15[i],2)+pow(err_dz3f_10_15[i]/y_dz3f_10_15[i],2));
    y_dz3f_10_15[i]/=y_dz10f_10_15[i];
  }
  TGraphErrors * r_dz3f_0_5 = new TGraphErrors(dz10f_0_5->GetN(),x_dz3f_0_5,y_dz3f_0_5,0,err_dz3f_0_5);
  TGraphErrors * r_dz3f_5_10 = new TGraphErrors(dz10f_0_5->GetN(),x_dz3f_5_10,y_dz3f_5_10,0,err_dz3f_5_10);
  TGraphErrors * r_dz3f_10_15 = new TGraphErrors(dz10f_0_5->GetN(),x_dz3f_10_15,y_dz3f_10_15,0,err_dz3f_10_15);
  r_dz3f_0_5->SetMarkerStyle(24);
  r_dz3f_0_5->SetMarkerColor(kBlue);

  r_dz3f_5_10->SetMarkerStyle(25);
  r_dz3f_5_10->SetMarkerColor(kRed);

  r_dz3f_10_15->SetMarkerStyle(30);
  r_dz3f_10_15->SetMarkerColor(kMagenta+1);


  TCanvas * can2 = new TCanvas("FakeCorrectionRatios","FakeCorrectionRatios",1300,750);
  can2->cd();
  TPad * pad3 = new TPad("pad3","pad3",0.0,0.0,xpad1,1.0);
  pad3->Draw();
  pad3->SetLeftMargin(lmargin1/xpad1);
  pad3->SetRightMargin(rmargin1/xpad1);
  pad3->cd();
  TH1D * hR = new TH1D("Ratio","Ratio",100,0,1.2);
  hR->SetMinimum(0.8);
  hR->SetMaximum(1.2);
  hR->SetXTitle("p_{T}");
  hR->SetYTitle("{#Delta z=6; #chi^{2}=20}/{#Delta z=10; #chi^{2}=40}");
  hR->SetStats(kFALSE);

  hR->Draw();
  gPad->SetGrid(1,1);
  r_dz3_0_5->Draw("p");
  r_dz3_5_10->Draw("p");
  r_dz3_10_15->Draw("p");
  TLegend * leg2 = new TLegend(0.24,0.6,0.44,0.9);
  leg2->AddEntry(r_dz3_0_5,   Form("0-5%c (3/20)",'%'),"p");
  leg2->AddEntry(r_dz3_5_10,  Form("5-10%c (3/20)",'%'),"p");
  leg2->AddEntry(r_dz3_10_15, Form("10-15%c (3/20)",'%'),"p");
  leg2->Draw();
  TLatex * text2 = new TLatex(0.7,0.012,"Uncorrected");
  text2->Draw();

  can2->cd();
  TPad * pad4 = new TPad("pad4","pad4",xpad2,0.,1.,1.0);
  pad4->Draw();
  pad4->SetLeftMargin(lmargin2/xpad2);
  pad4->SetRightMargin(rmargin2/xpad2);
  pad4->cd();
  hR->Draw();
  gPad->SetGrid(1,1);
  r_dz3f_0_5->Draw("p");
  r_dz3f_5_10->Draw("p");
  r_dz3f_10_15->Draw("p");
  can2->Print("FakeCorrectionRatio.png","png");


}
TGraphErrors * InputType1(TString fileName, Int_t &Num){
  cout<<"Open "<<fileName.Data()<<endl;
  FILE * fin = fopen(fileName.Data(),"r");
  Double_t x[100];
  Double_t y[100];
  Double_t dy[100];
  Num = 0;
  if(fin==NULL) {
    cout<<"File not found: "<<fileName.Data()<<endl;
    return 0;
  }
  char buf[80];
  fgets(buf,80,fin);
  while(fgets(buf,80,fin)!=NULL) {
    if(buf[0]=='_') continue;
    if(sscanf(buf,"%lf\t%lf\t%lf",&x[Num],&y[Num],&dy[Num]) == 3) {
      ++Num;
    } else {
      cout<<"Failed to interpret: "<<buf<<endl;
    };

  }
  if(Num==0) return 0;
  TGraphErrors * g = new TGraphErrors(Num,x,y,0,dy);
  g->SetName(fileName.Data());
  g->SetTitle(fileName.Data());
  return g;
}
