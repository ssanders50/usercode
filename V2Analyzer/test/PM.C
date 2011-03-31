void * PM(){
  TFile * tf = new TFile("EPSpectra.root");
  TGraphErrors * gm = (TGraphErrors *) tf->Get("-2.4_-2.0/dNdPt_40_50");
  TGraphErrors * gp = (TGraphErrors *) tf->Get("02.0_02.4/dNdPt_40_50");
  Double_t *pt=gm->GetX();
  Double_t *ym=gm->GetY();
  Double_t *yem=gm->GetEY();
  Double_t *yp=gp->GetY();
  Double_t *yep=gp->GetEY();
  Double_t xr[40];
  Double_t r[40];
  Double_t re[40];
  cout<<gm->GetN()<<endl;
  int n = 0;
  for(int i = 0; i<gm->GetN(); i++) {
    cout<<i<<" "<<ym[i]<<" "<<yp[i]<<" "<<yem[i]<<" "<<yep[i]<<endl;
    if(ym[i]>0) { 
      r[n] = yp[i]/ym[i];
      re[n] = sqrt(pow(yem[i]/ym[i],2)+pow(yep[i]/yp[i],2));
      xr[n++]=pt[i];
    }
  }
  TGraphErrors * gr = new TGraphErrors(n,xr,r,0,re);
  TH1D * hr = new TH1D("ratio","p/m",100,0,6);
  hr->SetXTitle("p_{T} (GeV/c)");
  hr->SetYTitle("Ratio(plus/minus)");
  hr->SetStats(kFALSE);
  hr->SetMaximum(1.4);
  hr->Draw();
  gPad->SetGrid(1,1);
  gr->SetMarkerStyle(22);
  gr->Draw("p");
  c1->Print("ratio.png","png");
}
