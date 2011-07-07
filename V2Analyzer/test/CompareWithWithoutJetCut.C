#include <iomanip>
#include <fstream>
TString TagLabel;
TGraphErrors * InputType1(TString fileName, Int_t & Num);
TGraphErrors * InputType1(int  minCent, int maxCent, double eta, TString tag, TString EP1, TString EP2,Int_t & Num) {
  TString name = "";
  if(EP2.Length()>0) {
    Form("results/v2pt_EP_%d-%d_%02d-%s_%s_%s.txt",minCent,maxCent,(int)(10*eta),tag.Data(),EP1.Data(),EP2.Data());
    TagLabel = Form("%s_%s_%s",tag.Data(),EP1.Data(),EP2.Data());
  } else {
    name = Form("results/v2pt_EP_%d-%d_%02d_%s_%s.txt",minCent,maxCent,(int)(10*eta),tag.Data(),EP1.Data());
    TagLabel = Form("%s_%s",tag.Data(),EP1.Data());
  }
  return InputType1(name,Num);
}
void CompareWithWithoutJetCut(){
  //                   0-5   5-10  10-15    15-20    20-25    25-30     30-35   35-40    40-50  50-60   60-70     70-80 80-90 90-100
  int mincent[14] ={     0,     5,    10,      15,      20,      25,      30,     35,       40,    50,     60,       70,   80,    90};
  int maxcent[14] ={     5,    10,    15,      20,      25,      30,      35,     40,       50,    60,     70,       80,   90,   100};
  int markers[14]= {    22,    21,    23,      24,      25,      30,      20,      3,       29,    28,     26,       27,    3,     5};
  int colors[14] = {kBlack, kBlue,  kRed, kCyan+2, kViolet, kSpring, kOrange, kRed+4, kAzure+9, kTeal, kRed-4,kYellow+2,    0,     0};


  TString intTag = "";
  TString tag = "v3v4.txt";


  Double_t lmargin1 = 0.1;
  Double_t rmargin1 = 0.01;
  Double_t lmargin2 = 0.1;
  Double_t rmargin2 = 0.05;
  Double_t xaxislen = (1.-(rmargin1+lmargin1+rmargin2+lmargin2))/2.;
  Double_t xpad1 = lmargin1+xaxislen+rmargin1;
  Double_t xpad2 = 1-xpad1;
  TCanvas * can = new TCanvas("CompareWithWithoutJetCut","Compare v_{2} With and Without Jet Cut",1300,750);
  can->cd();
  TPad * pad1 = new TPad("pad1","pad1",0.0,0.0,xpad1,1.0);
  pad1->Draw();
  pad1->SetLeftMargin(lmargin1/xpad1);
  pad1->SetRightMargin(rmargin1/xpad1);
  pad1->cd();
  TH1D * h = new TH1D("h","",500,0,11.999);
  h->SetMaximum(0.3);
  h->SetMinimum(0.0);
  h->SetXTitle("p_{T}");
  h->SetYTitle("v_{2}{EP}");
  h->SetStats(kFALSE);
  h->Draw();
  TGraphErrors * g[12];
  TGraphErrors * sg[12];
  TGraphErrors * r[12];
  TGraphErrors * sr[12];
  //With jets
  TString tag = "Asym-dz10chi40";
  g[0] = InputType1(0,  5,1.0, tag,  "etHF", "",    13);
  g[1] = InputType1(5, 10,1.0, tag,  "etHF", "",    13);
  g[2] = InputType1(10,15,1.0, tag,  "etHF", "",    13);
  g[3] = InputType1(15,20,1.0, tag,  "etHF", "",    13);
  g[4] = InputType1(20,25,1.0, tag,  "etHF", "",    13);
  g[5] = InputType1(25,30,1.0, tag,  "etHF", "",    13);
  g[6] = InputType1(30,35,1.0, tag,  "etHF", "",    13);
  g[7] = InputType1(35,40,1.0, tag,  "etHF", "",    13);
  g[8] = InputType1(40,50,1.0, tag,  "etHF", "",    13);
  g[9] = InputType1(50,60,1.0, tag,  "etHF", "",    13);

  Double_t xx[10][20];
  Double_t yy[10][20];
  Double_t dyy[10][20];
 
  for(Int_t i = 0; i<10; i++) {
    Double_t * x =  g[i]->GetX();
    Double_t * y =  g[i]->GetY();
    Double_t * dy = g[i]->GetEY();
    
    for(Int_t j=0; j<g[i]->GetN(); j++) {
      xx[i][j] = x[j];
      yy[i][j] = y[j];
      dyy[i][j] = dy[j];
    }
  }
 
  //Without jets
  TString tag = "Asym-dz10chi40jetCut1";
  sg[0] = InputType1(0,  5,1.0, tag,  "etHF", "",    13);
  sg[1] = InputType1(5, 10,1.0, tag,  "etHF", "",    13);
  sg[2] = InputType1(10,15,1.0, tag,  "etHF", "",    13);
  sg[3] = InputType1(15,20,1.0, tag,  "etHF", "",    13);
  sg[4] = InputType1(20,25,1.0, tag,  "etHF", "",    13);
  sg[5] = InputType1(25,30,1.0, tag,  "etHF", "",    13);
  sg[6] = InputType1(30,35,1.0, tag,  "etHF", "",    13);
  sg[7] = InputType1(35,40,1.0, tag,  "etHF", "",    13);
  sg[8] = InputType1(40,50,1.0, tag,  "etHF", "",    13);
  sg[9] = InputType1(50,60,1.0, tag,  "etHF", "",    13);


  for(Int_t i = 0; i<10; i++) {
    g[i]->SetMarkerStyle(20);
    g[i]->SetMarkerColor(colors[i]);
    g[i]->SetLineColor(colors[i]);
    g[i]->SetMarkerSize(1.2);
    g[i]->Draw("P");
  }

  for(Int_t i = 0; i<10; i++) {
    sg[i]->SetMarkerStyle(25);
    sg[i]->SetMarkerColor(colors[i]);
    sg[i]->SetLineColor(colors[i]);
    sg[i]->SetMarkerSize(1.2);
    sg[i]->Draw("P");
  }

  TLegend * leg = new TLegend(0.82,0.66,0.94,0.93,"");
  leg->SetFillColor(0);
  for(Int_t i = 9; i>=0; i--) {
    leg->AddEntry(g[i],  Form( "%d-%d%c",mincent[i],maxcent[i],'%'),"lp");
  }
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0.0);
  leg->Draw();


  can->cd();
  TPad * pad2 = new TPad("pad2","pad2",xpad1,0.0,1.0,1.0);
  pad2->Draw();
  pad2->SetLeftMargin(lmargin2/xpad2);
  pad2->SetRightMargin(rmargin2/xpad2);
  pad2->cd();
  
  gPad->SetGrid(1,1);
  TH1D * hr = new TH1D("hr","",500,0,11.999);
  hr->SetMaximum(1.2);
  hr->SetMinimum(0.0);
  hr->SetYTitle("Ratio (Jet Cut/All)");
  hr->SetXTitle("p_{T}");
  hr->Draw();
  
  for(int i = 0; i<10; i++) {
    Double_t * x = sg[i]->GetX();
    Double_t * y = sg[i]->GetY();
    Double_t * dy = sg[i]->GetEY();
    for(int j = 0; j<sg[i]->GetN(); j++) {
     yy[i][j] = y[j]/yy[i][j];
     double e1 = dyy[i][j]/yy[i][j];
     double e2 = dy[j]/y[j];
     double sq = TMath::Sqrt(e1*e1+e2*e2);
     dyy[i][j] =yy[i][j]*sq;
   }
   r[i] = new TGraphErrors(sg[i]->GetN(),xx[i],yy[i],0,dyy[i]);
   r[i]->SetMarkerStyle(20);
   r[i]->SetMarkerColor(colors[i]);
   r[i]->SetLineColor(colors[i]);
   r[i]->SetLineWidth(2);
   r[i]->SetMarkerSize(1.4);
   r[i]->Draw("P");
   }
  can->Print("WithWithoutJetCut.png","png");

}
TGraphErrors * InputType1(TString fileName, Int_t &Num){
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
      //      cout<<"Failed to interpret: "<<buf<<endl;
    };

  }
  if(Num==0) return 0;
  TGraphErrors * g = new TGraphErrors(Num,x,y,0,dy);
  g->SetName(fileName.Data());
  g->SetTitle(fileName.Data());
  return g;
}
