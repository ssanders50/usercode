
#include <iomanip>
#include <fstream>

void makeNormalizedEtaDist(){
  TString tag = "Flow_Skim_Run2010-v5_dz5Flat_-10to10.txt";
  static const int neta = 6;
  TFile * tf[neta];
  //double etabins[] = {-2.4,-2.0,-1.6,-0.8,0,0.8,1.6,2.0,2.4};
  double etabins[]={ 0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
  //                   0-5   5-10  10-15    15-20    20-25    25-30     30-35   35-40    40-50  50-60   60-70     70-80 80-90 90-100
  int mincent[14] ={     0,     5,    10,      15,      20,      25,      30,     35,       40,    50,     60,       70,   80,    90};
  int maxcent[14] ={     5,    10,    15,      20,      25,      30,      35,     40,       50,    60,     70,       80,   90,   100};
  int markers[14]= {    22,    21,    23,      24,      25,      30,      20,      3,       29,    28,     26,       27,    3,     5};
  int colors[14] = {kBlack, kBlue,  kRed, kCyan+2, kViolet, kSpring, kOrange, kRed+4, kAzure+9, kTeal, kRed-4,kYellow+2,    0,     0};
  double y_0_5[neta];
  double y_5_10[neta];
  double y_10_15[neta];
  double y_15_20[neta];
  double y_20_25[neta];
  double y_25_30[neta];
  double y_30_35[neta];
  double y_35_40[neta];
  double y_40_50[neta];
  double y_50_60[neta];
  double y_60_70[neta];
  double dy_0_5[neta];
  double dy_5_10[neta];
  double dy_10_15[neta];
  double dy_15_20[neta];
  double dy_20_25[neta];
  double dy_25_30[neta];
  double dy_30_35[neta];
  double dy_35_40[neta];
  double dy_40_50[neta];
  double dy_50_60[neta];
  double dy_60_70[neta];
  TF1 * fit[12];
  double x[neta] = {0.2,0.6,1.0,1.4,1.8,2.2};
  FILE * fin ;
  char buf[80];
  TString intTag = "";
  TCanvas * can = new TCanvas("IntV2Eta","IntV2 vs. Eta",1100,750);
  can->Divide(2);
  can->cd(1);
  gPad->SetLeftMargin(0.21);
  for(int i = 0; i< neta; i++ ) {
    TString fname = Form("results/intv2%s_%03d_%03d_%03d_%03d_%s",intTag.Data(),
			 (int)(-10*etabins[i+1]),(int)(-10*etabins[i]),(int)(10*etabins[i]),(int)(10*etabins[i+1]),tag.Data());
    cout<<fname.Data()<<endl;
    fin = fopen(fname.Data(),"r");
    double tmp;
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_0_5[i],&dy_0_5[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_5_10[i],&dy_5_10[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_10_15[i],&dy_10_15[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_15_20[i],&dy_15_20[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_20_25[i],&dy_20_25[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_25_30[i],&dy_25_30[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_30_35[i],&dy_30_35[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_35_40[i],&dy_35_40[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_40_50[i],&dy_40_50[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_50_60[i],&dy_50_60[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_60_70[i],&dy_60_70[i]);
    fclose(fin);
  }
  TH1D * h = new TH1D("h","",500,0,3.5);
  h->SetMaximum(0.12);
  h->SetMinimum(0.0);
  h->SetXTitle("#eta");
  h->SetYTitle("v_{2}{EP}");
  h->SetStats(kFALSE);
  h->GetYaxis()->SetTitleOffset(1.7);
  h->Draw();
  gPad->SetGrid(1,1);
  std::ofstream file;
  TGraphErrors * g[11];
  g[0] = new TGraphErrors(neta,x,y_0_5,0,  dy_0_5);
  g[1] = new TGraphErrors(neta,x,y_5_10,0, dy_5_10);
  g[2] = new TGraphErrors(neta,x,y_10_15,0,dy_10_15);
  g[3] = new TGraphErrors(neta,x,y_15_20,0,dy_15_20);
  g[4] = new TGraphErrors(neta,x,y_20_25,0,dy_20_25);
  g[5] = new TGraphErrors(neta,x,y_25_30,0,dy_25_30);
  g[6] = new TGraphErrors(neta,x,y_30_35,0,dy_30_35);
  g[7] = new TGraphErrors(neta,x,y_35_40,0,dy_35_40);
  g[8] = new TGraphErrors(neta,x,y_40_50,0,dy_40_50);
  g[9] = new TGraphErrors(neta,x,y_50_60,0,dy_50_60);
  g[10] = new TGraphErrors(neta,x,y_60_70,0,dy_60_70);
  TLegend * leg = new TLegend(0.75,0.5,0.92,0.89,"Centrality");
  Double_t yslope[11];
  Double_t dyslope[11];
  Double_t xslope[11];
  for(Int_t i = 0; i<11; i++) {
    TString fileName = Form("results/intV2_EP_%d-%d.txt",mincent[i],maxcent[i]);
    file.open(fileName.Data());
    g[i]->SetMarkerStyle(markers[i]);
    g[i]->SetMarkerColor(colors[i]);
    g[i]->SetLineColor(colors[i]);
    g[i]->Draw("P");
    fit[i] = new TF1(Form("fit_%d",i),"pol1",0.1,2.4);
    fit[i]->SetLineColor(colors[i]);
    fit[i]->SetLineStyle(2);
    g[i]->Fit(fit[i],"QR");
    xslope[i]=(mincent[i]+maxcent[i])/2.;
    yslope[i]=fabs(fit[i]->GetParameter(1));
    dyslope[i]=fit[i]->GetParError(1);
    cout<<xslope[i]<<" "<<yslope[i]<<" "<<dyslope[i]<<endl;
    Double_t * xx = g[i]->GetX();
    Double_t * yy = g[i]->GetY();
    Double_t * eyy = g[i]->GetEY();
    leg->AddEntry(g[i],  Form( "%d-%d",mincent[i],maxcent[i]),"lp");
    for(int j = 0; j< neta; j++) {
      file<<setprecision(3)<<xx[j]<<"\t"<<setprecision(3)<<yy[j]<<"\t"<<setprecision(3)<<eyy[j]<<"\t"<<endl;
    }
    file.close();
  }
  leg->SetTextSize(0.03);
  leg->Draw();
  //TLatex * text = new TLatex(-2.5,0.15,"Integration Range: 0.3 - 3 GeV/c");
  //TLatex * teff = new TLatex(-2.5,0.13,Form("Efficiency Tables: %s",intTag.Data()));
  //text->Draw();
  //teff->Draw();
  can->cd(2);
  gPad->SetLeftMargin(0.21);
  TGraphErrors * gslope = new TGraphErrors(11,xslope,yslope,0,dyslope);
  TH1D * hslope = new TH1D("hslope","",200,0,100);
  hslope->SetStats(kFALSE);
  hslope->SetMaximum(0.02);
  hslope->SetXTitle("Centrality");
  hslope->SetYTitle("|dv_{2}/d#eta|");
  hslope->GetYaxis()->SetTitleOffset(2.0);
  hslope->Draw();
  gslope->SetMarkerStyle(21);
  gslope->SetMarkerColor(kBlue);
  gslope->SetLineColor(kBlue);
  gslope->Draw("p");
  can->Print(Form("v2Eta_%s.png",intTag.Data()));
}
