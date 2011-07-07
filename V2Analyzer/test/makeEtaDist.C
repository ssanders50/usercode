#include <iomanip>
#include <fstream>

void makeEtaDist(){
  TString tag = "_hiGoodMergedTracks_03_2PtCut_PtWeight.txt";
  TFile * tf[8];
  double etabins[] = {-2.4,-2.0,-1.6,-0.8,0,0.8,1.6,2.0,2.4};
  //                0-5   5-10 10-15 15-20   20-30   30-40   40-50   50-60  60-70    70-80 80-90 90-100
  int mincent[12] ={0,     5,    10,  15,     20,     30,     40,     50,    60,      70,   80,    90       };
  int maxcent[12] ={5,    10,    15,  20,     30,     40,     50,     60,    70,      80,   90,   100       };
  int markers[12]= {22,    21,   23,  24,     25,     20,     29,     28,    26,      27,   3,     5        };
  int colors[12] = {kBlack,kBlue,kRed,kCyan+2,kViolet,kSpring,kOrange,kRed+4,kAzure+9,kTeal,kRed-4,kYellow+2};
  double y_0_5[8];
  double y_5_10[8];
  double y_10_15[8];
  double y_15_20[8];
  double y_20_30[8];
  double y_30_40[8];
  double y_40_50[8];
  double y_50_60[8];
  double y_60_70[8];
  double dy_0_5[8];
  double dy_5_10[8];
  double dy_10_15[8];
  double dy_15_20[8];
  double dy_20_30[8];
  double dy_30_40[8];
  double dy_40_50[8];
  double dy_50_60[8];
  double dy_60_70[8];
  double x[8] = {-2.2,-1.8,-1.2,-0.4,0.4,1.2,1.8,2.2};
  FILE * fin ;
  char buf[80];
  for(int i = 0; i< 8; i++ ) {
    TString fname = Form("results/intv2_%03d_%03d_%s",(int)(10*etabins[i]),(int)(10*etabins[i+1]),tag.Data());
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
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_20_30[i],&dy_20_30[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_30_40[i],&dy_30_40[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_40_50[i],&dy_40_50[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_50_60[i],&dy_50_60[i]);
    fgets(buf,80,fin);
    sscanf(buf,"%lf\t%lf\t%lf\n",&tmp,&y_60_70[i],&dy_60_70[i]);
    fclose(fin);
  }
  TH1D * h = new TH1D("h","",500,-3,3);
  h->SetMaximum(0.25);
  h->SetMinimum(0.0);
  h->SetXTitle("#eta");
  h->SetYTitle("v_{2}{RP}");
  h->SetStats(kFALSE);
  h->Draw();
  gPad->SetGrid(1,1);
  std::ofstream file;
  TGraphErrors * g[9];
  g[0] = new TGraphErrors(8,x,y_0_5,0,  dy_0_5);
  g[1] = new TGraphErrors(8,x,y_5_10,0, dy_5_10);
  g[2] = new TGraphErrors(8,x,y_10_15,0,dy_10_15);
  g[3] = new TGraphErrors(8,x,y_15_20,0,dy_15_20);
  g[4] = new TGraphErrors(8,x,y_20_30,0,dy_20_30);
  g[5] = new TGraphErrors(8,x,y_30_40,0,dy_30_40);
  g[6] = new TGraphErrors(8,x,y_40_50,0,dy_40_50);
  g[7] = new TGraphErrors(8,x,y_50_60,0,dy_50_60);
  g[8] = new TGraphErrors(8,x,y_60_70,0,dy_60_70);
  TLegend * leg = new TLegend(0.75,0.67,0.92,0.92,"Centrality");
  for(Int_t i = 0; i<9; i++) {
    TString fileName = Form("results/intV2_EP_%d-%d.txt",mincent[i],maxcent[i]);
    file.open(fileName.Data());
    g[i]->SetMarkerStyle(markers[i]);
    g[i]->SetMarkerColor(colors[i]);
    g[i]->SetLineColor(colors[i]);
    g[i]->Draw("P");
    Double_t * xx = g[i]->GetX();
    Double_t * yy = g[i]->GetY();
    Double_t * eyy = g[i]->GetEY();
    leg->AddEntry(g[i],  Form( "%d-%d",mincent[i],maxcent[i]),"lp");
    for(int j = 0; j< 8; j++) {
      file<<setprecision(3)<<xx[j]<<"\t"<<setprecision(3)<<yy[j]<<"\t"<<setprecision(3)<<eyy[j]<<"\t"<<endl;
    }
    file.close();
  }
  leg->Draw();
  TLatex * text = new TLatex(-2.5,0.2,"Integration Range: 0.3 - 3 GeV/c");
  text->Draw();
  c1->Print("~/public_html/V2Spectra/v2Eta.png");
}
