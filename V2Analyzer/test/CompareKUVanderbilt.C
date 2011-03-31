TFile * tf;
static const Int_t nPtBins = 15;
static const double ptbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0};
void CompareKUVanderbilt() {
  tf = new TFile("EPSpectra.root");
  TGraphErrors * ku5_10 = tf->Get("-0.8_00.0_00.0_00.8/gd2NdEtadPt_5_10");
  TH1D * hpt = new TH1D("hpt","hpt",nPtBins,ptbins);
  TH1D * hist = (TH1D *)ku5_10->GetHistogram();
  hist->SetMinimum(0.001);
  hist->SetMaximum(10000.);
  ku5_10->GetHistogram()->Draw();
  gPad->SetLogy();
  ku5_10->Draw("p");
  double * xx = ku5_10->GetX();
  double * yy = ku5_10->GetY();
  double sumKU = 0;
  for(int i = 0; i<ku5_10->GetN(); i++) {
    double wKU = hpt->GetBinWidth(i+1);
    if(xx[i]>ptbins[1]&&xx[i]<ptbins[11]) {
      sumKU+=yy[i]*wKU;
      cout<<xx[i]<<" "<<yy[i]<<" "<<wKU<<" -- "<<sumKU<<endl;
    }
  }
  cout<<" "<<endl;
  FILE * fin;
  fin = fopen("dNdpt_LYZ/dNdpt_LYZ_0510.txt","r");
  double ptV[14];
  double dndptV[14];
  double f1,f2,f3,f4,f5;
  double sumV = 0;
  for(int i = 0; i<14; i++) {
    fscanf(fin,"%lf %lf %lf %lf %lf",&f1,&f2,&f3,&f4,&f5);
    ptV[i]=f1;
    dndptV[i]=f5;
    double wV = hpt->GetBinWidth(i+2);
    if(f1>ptbins[1]&&f1<ptbins[11]) {
      sumV+=dndptV[i]*wV;
      cout<<ptV[i]<<" "<<dndptV[i]<<" "<<wV<<" -- "<<sumV<<endl;
    } 
 }
  TGraph * Van = new TGraph(14,ptV,dndptV);
  Van->SetMarkerStyle(21);
  Van->SetMarkerColor(kRed);
  Van->Draw("p");

}
