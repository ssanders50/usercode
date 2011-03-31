// -*- C++ -*-
//
// Package:    V2Analyzer
// Class:      V2Analyzer
// 
/* class V2Analyzer V2Analyzer.cc HiEvtPlaneFlatten/V2Analyzer/src/V2Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Stephen Sanders
//         Created:  Wed Jun 23 12:27:13 EDT 2010
// $Id: V2Analyzer.cc,v 1.3 2010/07/26 23:11:00 ssanders Exp $
//
//

#define TRACKCOLLECTION 1
//#define RECOCHARGEDCANDIDATECOLLECTION 1



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/HeavyIon.h"
#include "HepMC/SimpleVector.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"

#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/UpdatablePSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "CondFormats/DataRecord/interface/RPFlatParamsRcd.h"
#include "CondFormats/RPFlatParams/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include <time.h>
#include <cstdlib>

using namespace std;
#include <vector>
using std::vector;
using std::rand;
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
static const double pi = 3.14159265358979312;
static const double pi2 = 1.57079632679489656;

static const Int_t nCentBins = 14;
static const Int_t nPtBins = 15;
static const Int_t nEtBins = 15;
//static const Int_t nEtaBins = 22;
static const Int_t nEtaBins = 50;
static const double centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,90,100};
static const double ptbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0,16.};
static const double etbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0,16.};
//static const double etabins[]={-5,-4.5,-4,-3.5,-3,-2.4,-2,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2.0,2.4,3,3.5,4,4.5,5};
static const double etabins[]={-5.0, -4.8, -4.6, -4.4, -4.2, -4.0, -3.8, -3.6, -3.4, -3.2,
			       -3.0, -2.8, -2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2,
			       -1.0, -0.8, -0.6, -0.4, -0.2,  0.0,  0.2,  0.4,  0.6,  0.8,
                                1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,  2.8,
			        3.0,  3.2,  3.4,  3.6,  3.8,  4.0,  4.2,  4.4,  4.6,  4.8,
			        5.0};

//
// class declaration
//

class V2Analyzer : public edm::EDAnalyzer {
public:
  explicit V2Analyzer(const edm::ParameterSet&);
  ~V2Analyzer();
 
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  double dzerr_;
  double chi2_; 
  class v2Generator {
  private:
  Double_t bounds2(Double_t ang) {
    if(ang<-pi2) ang+=pi;
    if(ang>pi2)  ang-=pi;
    return ang;
  }
  public:
    explicit v2Generator(TFileDirectory  dir, string label, double gap){
      label_ = label;
      gap_ = fabs(gap);
      if(gap_ == 0) {
	subAuto = false;
      } else {
	subAuto = true;
	static const int maxAutoEtaBins = 4000;
	int numAutoEtaBins = (etabins[nEtaBins]-etabins[0])/gap_;
	if(numAutoEtaBins > maxAutoEtaBins) {
	  gap_ = (etabins[nEtaBins]-etabins[0])/maxAutoEtaBins;
	  numAutoEtaBins = maxAutoEtaBins;
	}
	autoSin = dir.make<TH1D>(Form("autoSin_%s",label.data()),
				 Form("autoSin_%s",label.data()),4*numAutoEtaBins,etabins[0],etabins[nEtaBins]);
	autoSin->Sumw2();
	autoCos = dir.make<TH1D>(Form("autoCos_%s",label.data()),
				 Form("autoCos_%s",label.data()),4*numAutoEtaBins,etabins[0],etabins[nEtaBins]);
	autoCos->Sumw2();
 	autoCnt = dir.make<TH1D>(Form("autoCnt_%s",label.data()),
				 Form("autoCnt_%s",label.data()),4*numAutoEtaBins,etabins[0],etabins[nEtaBins]);
	autoCnt->Sumw2();
	autoPsi = dir.make<TH1D>(Form("autoPsi_%s",label.data()),Form("autoPsi_%s",label.data()),200,-4,4);
	autoPsi->Sumw2();
      }
      eta_ = dir.make<TH1F>(Form("eta_%s",label.data()),Form("eta_%s",label.data()),nEtaBins,etabins);
      for(int i = 0; i<nEtaBins; i++) {
	cos_[i] = dir.make<TH2D>(Form("cos_%s_%d",label.data(),i),Form("cos_%s_%d",label.data(),i),nCentBins,centbins,nPtBins,ptbins);
	cos_[i]->SetXTitle("cent");
	cos_[i]->SetYTitle("p_{T}");
	cos_[i]->SetOption("colz");
	cos_[i]->Sumw2();
	psi_[i] = dir.make<TH2D>(Form("psi_%s_%d",label.data(),i),Form("psi_%s_%d",label.data(),i),nCentBins,centbins,50,-2,2);
	psi_[i]->SetXTitle("cent");
	psi_[i]->SetYTitle("#psi");
	psi_[i]->SetOption("colz");
	psi_[i]->Sumw2();
	genpsi_[i] = dir.make<TH2D>(Form("genpsi_%s_%d",label.data(),i),Form("genpsi_%s_%d",label.data(),i),nCentBins,centbins,50,-2,2);
	genpsi_[i]->SetXTitle("cent");
	genpsi_[i]->SetYTitle("#psi");
	genpsi_[i]->SetOption("colz");
	genpsi_[i]->Sumw2();
	cnt_[i] = dir.make<TH2D>(Form("cnt_%s_%d",label.data(),i),Form("cnt_%s_%d",label.data(),i),nCentBins,centbins,nPtBins,ptbins);
	cnt_[i]->SetXTitle("cent");
	cnt_[i]->SetYTitle("p_{T}");
	cnt_[i]->SetOption("colz");
	cnt_[i]->Sumw2();
      } 
      for(int i = 0; i<nEtaBins; i++) {
	gencos_[i] = dir.make<TH2D>(Form("gencos_%s_%d",label.data(),i),Form("gencos_%s_%d",label.data(),i),nCentBins,centbins,nPtBins,ptbins);
	gencos_[i]->SetXTitle("cent");
	gencos_[i]->SetYTitle("p_{T}");
	gencos_[i]->SetOption("colz");
	gencos_[i]->Sumw2();
	gencnt_[i] = dir.make<TH2D>(Form("gencnt_%s_%d",label.data(),i),Form("gencnt_%s_%d",label.data(),i),nCentBins,centbins,nPtBins,ptbins);
	gencnt_[i]->SetXTitle("cent");
	gencnt_[i]->SetYTitle("p_{T}");
	gencnt_[i]->SetOption("colz");
	gencnt_[i]->Sumw2();
      } 
      
      c12_ = dir.make<TH1D>(Form("c12_%s",label.data()),Form("c12_%s",label.data()),nCentBins,centbins);
      c12_->Sumw2();
      c13_ = dir.make<TH1D>(Form("c13_%s",label.data()),Form("c13_%s",label.data()),nCentBins,centbins);
      c13_->Sumw2();
      c23_ = dir.make<TH1D>(Form("c23_%s",label.data()),Form("c23_%s",label.data()),nCentBins,centbins);
      c23_->Sumw2();
      ccnt12_ = dir.make<TH1D>(Form("ccnt12_%s",label.data()),Form("ccnt12_%s",label.data()),nCentBins,centbins);
      ccnt12_->Sumw2();
      ccnt13_ = dir.make<TH1D>(Form("ccnt13_%s",label.data()),Form("ccnt13_%s",label.data()),nCentBins,centbins);
      ccnt13_->Sumw2();
      ccnt23_ = dir.make<TH1D>(Form("ccnt23_%s",label.data()),Form("ccnt23_%s",label.data()),nCentBins,centbins);
      ccnt23_->Sumw2();
      genres = dir.make<TH1D>(Form("genres_%s",label.data()),Form("genres_%s",label.data()),nCentBins,centbins);
      genres->Sumw2();
      genrescnt = dir.make<TH1D>(Form("genrescnt_%s",label.data()),Form("genrescnt_%s",label.data()),nCentBins,centbins);
      genrescnt->Sumw2();
      hphi = dir.make<TH1D>(Form("phi_%s",label.data()),Form("phi_%s",label.data()),100,-4,4);
      hPsi = dir.make<TH1D>(Form("Psi_%s",label.data()),Form("Psi_%s",label.data()),100,-4,4);
      hPsiGen = dir.make<TH1D>(Form("PsiGen_%s",label.data()),Form("PsiGen_%s",label.data()),100,-4,4);
   

    }
    ~v2Generator() ;
    
 
    void AddParticle(double phi, double Psi, double cent, double eta, double pt) {
      int ietabin = eta_->FindBin(eta)-1;
      if(ietabin<0 || ietabin>=nEtaBins) return;
      if(Psi<-4) return;
      cos_[ietabin]->Fill(cent,pt,TMath::Cos(2*(phi-Psi)));
      psi_[ietabin]->Fill(cent,bounds2(phi-Psi));
      cnt_[ietabin]->Fill(cent,pt);
      hphi->Fill(phi);
    }
    void AddGenParticle(double phi, double Psi, double cent, double eta, double pt) {
      int ietabin = eta_->FindBin(eta)-1;
      if(ietabin<0 || ietabin>=nEtaBins) return;
      if(Psi<-4) return;
      gencos_[ietabin]->Fill(cent,pt,TMath::Cos(2*(phi-Psi)));
      genpsi_[ietabin]->Fill(cent,bounds2(phi-Psi));
      gencnt_[ietabin]->Fill(cent,pt);
    }
    void SetAutocorrelation(double phi, double eta, double w) {
      if(!subAuto) return;
      if(w<0.1) return;
      autoSin->Fill(eta, w*sin(2. * phi));
      autoCos->Fill(eta, w*cos(2. * phi));
      autoCnt->Fill(eta);
    }
    void ResetAutocorrelation() {
      if(! subAuto) return;
      autoSin->Reset();
      autoCos->Reset();
      autoCnt->Reset();
    }
    void AddToResCor(double phiA, double phiB, double phiC, double cent) {
      if(phiA>-4 && phiB>-4) {
	c12_->Fill(cent,  TMath::Cos(2*(phiA-phiB)));
	ccnt12_->Fill(cent);
      }
      if(phiA>-4 && phiC>-4) {
	c13_->Fill(cent, TMath::Cos(2*(phiA-phiC)));
	ccnt13_->Fill(cent);
      }
      if(phiB>-4 && phiC>-4) {
	c23_->Fill(cent, TMath::Cos(2*(phiB-phiC)));
	ccnt23_->Fill(cent);
      }
      if(phiA>-4) hPsi->Fill(phiA);
    }
    
    void AddToGenRes(double phi, double genphi, double cent) {
      if(phi<-5) return;
      genres->Fill(cent,TMath::Cos(2*(phi-genphi)));
      genrescnt->Fill(cent);
      hPsiGen->Fill(genphi);
    }

    double GetAutoCorrectedPsi(double eta, double psi, double fSin, double fCos) {
      if(!subAuto) return psi;
      int lbinMin = autoCnt->FindBin(eta-gap_);
      int lbinMax = autoCnt->FindBin(eta+gap_);
      int multLost = 0;
      for(int i = lbinMin; i<=lbinMax; i++) {
	fSin-=autoSin->GetBinContent(i);
	fCos-=autoCos->GetBinContent(i);
	multLost+=autoCnt->GetBinContent(i);
      }
      psi = 0.5*atan2(fSin,fCos);
      autoPsi->Fill(psi);
      return psi;
    }
    TH2D * GetCos(int indx){return cos_[indx];}
    TH2D * GetCnt(int indx){return cnt_[indx];}
  private:
    TDirectory * saveDir;
    string label_;
    double gap_;
    TH2D * cos_[nEtaBins];
    TH2D * psi_[nEtaBins];
    TH2D * genpsi_[nEtaBins];
    TH2D * cnt_[nEtaBins];
    TH2D * gencos_[nEtaBins];
    TH2D * gencnt_[nEtaBins];
    int nEtaBins_; 
    TH1F * eta_;
    TH1D * autoSin;
    TH1D * autoCos;
    TH1D * autoCnt;
    TH1D * autoPsi; 
    TH1D * c12_;
    TH1D * c13_;
    TH1D * c23_;
    TH1D * ccnt12_;
    TH1D * ccnt13_;
    TH1D * ccnt23_;
    TH1D * genres;
    TH1D * genrescnt;
    TH1D * hphi;
    TH1D * hPsi;
    TH1D * hPsiGen;
    TH2D * hHiM;
    bool subAuto;
  };
  
  edm::Service<TFileService> fs;
  CentralityProvider * centrality_;
  //  const CentralityBins * cbins_;
  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TH2D * PhiEta_[nPtBins];
 TH2D * hTRACKp_etHFp;
  TH2D * hTRACKp_HFp;
  TH2D * hTRACKm_HFm;
  TH2D * hTRACKm_etHFm;
  TH2D * hTRACKp_TRACKm;
  TH2D * hetCaloP_etCaloM;
  TH2D * hetHFp_etHFm;
  TH2D * hetHFp_etHFmVsCent;
  TH2D * hMid_etHFpVsCent;
  TH2D * hPsi_GenPsi[NumEPNames];
  TH1D * hcent;
  TH1D * heta;
  TH1D * hMultByNpart;
  TH1D * hMultByNpartCnt;
  TH1D * hFull[NumEPNames];
  TH1D * hFullBin[NumEPNames][nCentBins];
  TH1D * hSub1Bin[NumEPNames][nCentBins];
  TH1D * hSub2Bin[NumEPNames][nCentBins];
  TH2D * hSub1Sub2[NumEPNames];
  TH1D * hGenRes[NumEPNames][nCentBins];
  TH1D * hSubRes[NumEPNames][nCentBins];
  TH1D * hMult1[NumEPNames];
  TH1D * hMult2[NumEPNames];
  TH1D * hMult[NumEPNames];
  TH1D * hMult1Cnt[NumEPNames];
  TH1D * hMult2Cnt[NumEPNames];
  TH1D * hMultCnt[NumEPNames];
  TH1D * hEta[NumEPNames];
  TH1D * hq[NumEPNames][nCentBins];
  TH1D * hNpartBin;
  TH1D * hNpartBinCnt;
  TH2D * hpt[nEtaBins];
  TH2D * hptCnt[nEtaBins];
  TH2D * het[nEtaBins];
  TH2D * hetCnt[nEtaBins];
  TH2D * hEmHad_EmEt[nEtaBins];
  TH2D * hPsiFlatVsPsiRaw[8];
  TH1D * hCentBinned;
  TH1D * hFlatDiffMean[NumEPNames][nCentBins];

  v2Generator * v2_Tracks[NumEPNames];
  v2Generator * v2_Calo[NumEPNames];
  Double_t bounds(Double_t ang) {
    if(ang<-pi) ang+=2.*pi;
    if(ang>pi)  ang-=2.*pi;
    return ang;
  }
  Double_t bounds2(Double_t ang) {
    if(ang<-pi2) ang+=pi;
    if(ang>pi2)  ang-=pi;
    return ang;
  }
  Bool_t InEtaRange(Int_t EP, Double_t eta) {
    if((eta>=EPEtaMin1[EP] && eta<EPEtaMax1[EP]) ||
       (EPEtaMin2[EP]!=EPEtaMax2[EP] && eta>=EPEtaMin2[EP] && eta<EPEtaMax2[EP])) {
      return kTRUE;
    } else {
      return kFALSE;
    }
  }
};

//
// constants, enums and typedefs
//
typedef std::vector<TrackingParticle>                   TrackingParticleCollection;
typedef TrackingParticleRefVector::iterator               tp_iterator;

//
// static data member definitions
//

//
// constructors and destructor
//
V2Analyzer::V2Analyzer(const edm::ParameterSet& iConfig)
  
{
  dzerr_ = iConfig.getUntrackedParameter<double>("dzerr_",14.);
  chi2_  = iConfig.getUntrackedParameter<double>("chi2_",80.);
  //now do what ever initialization is needed
  
  //  cbins_ = 0;
  centrality_ = 0;
  hcent = fs->make<TH1D>("cent","cent",200,-10,110);
  hCentBinned = fs->make<TH1D>("centBinned","centBinned",nCentBins,centbins);
  heta = fs->make<TH1D>("heta","heta",nEtaBins,etabins);
  hNpartBin = fs->make<TH1D>("NpartBin","NpartBin",nCentBins,centbins);
  hNpartBin->Sumw2();
  hNpartBinCnt = fs->make<TH1D>("NpartBinCnt","NpartBinCnt",nCentBins,centbins);
  hNpartBinCnt->Sumw2();
  TFileDirectory calodir = fs->mkdir("Calo");
  TFileDirectory specdir = fs->mkdir("Spectra");
  TFileDirectory rpcorrdir = fs->mkdir("RPCorrelations");
  TFileDirectory flatchkdir = fs->mkdir("FlatteningChecks");
  for(int i = 0; i< 8; i++) {
    hPsiFlatVsPsiRaw[i] = fs->make<TH2D>(Form("PhiFlatVsPhiRaw_%s",EPNames[i].data()),Form("PhiFlatVsPhiRaw_%s",EPNames[i].data()),50,-2,2,50,-2,2);
    hPsiFlatVsPsiRaw[i]->SetOption("colz");
    hPsiFlatVsPsiRaw[i]->SetXTitle("#Psi_{no flattening}");
    hPsiFlatVsPsiRaw[i]->SetYTitle("#Psi_{flattened}");
    hPsiFlatVsPsiRaw[i]->SetStats(kFALSE);
  }
  for(int i = 0; i< nEtaBins; i++) {
    hpt[i] = specdir.make<TH2D>(Form("pt_%d",i),Form("pt_%d",i),nPtBins,ptbins,nCentBins,centbins);
    hptCnt[i] = specdir.make<TH2D>(Form("ptCnt_%d",i),Form("ptCnt_%d",i),nPtBins,ptbins,nCentBins,centbins);
    hpt[i]->Sumw2();
    hptCnt[i]->Sumw2();
    hpt[i]->SetOption("colz");
    hpt[i]->SetXTitle("pt (GeV/c)");
    hpt[i]->SetYTitle("centrality");
    hpt[i]->SetTitle(Form("%5.2f #leq #eta < %5.2f",etabins[i],etabins[i+1]));
    hptCnt[i]->SetOption("colz");
    het[i] = specdir.make<TH2D>(Form("et_%d",i),Form("et_%d",i),nEtBins,etbins,nCentBins,centbins);
    hetCnt[i] = specdir.make<TH2D>(Form("etCnt_%d",i),Form("etCnt_%d",i),nEtBins,etbins,nCentBins,centbins);
    het[i]->Sumw2();
    het[i]->SetXTitle("et (GeV)");
    het[i]->SetYTitle("centrality");
    het[i]->SetTitle(Form("%5.2f #leq #eta < %5.2f",etabins[i],etabins[i+1]));
    hetCnt[i]->Sumw2();
    het[i]->SetOption("colz");
    hetCnt[i]->SetOption("colz");
    hEmHad_EmEt[i] = calodir.make<TH2D>(Form("EmHa_EmdEt_%d",i),Form("EmHad_EmEt_%d",i),60,-5,15,60,-5,15);
    hEmHad_EmEt[i]->SetTitle(Form("%5.2f #leq #eta < %5.2f",etabins[i],etabins[i+1]));
    hEmHad_EmEt[i]->SetXTitle("etEM+etHad (GeV)");
    hEmHad_EmEt[i]->SetYTitle("etEM (GeV)");
    hEmHad_EmEt[i]->SetOption("colz");
  }
  hMultByNpart = fs->make<TH1D>("MultByNpart","2*Mult/N_{part}",200,-10,110);
  hMultByNpartCnt = fs->make<TH1D>("MultByNpartCnt","MultByNpartCnt",200,-10,110);
  hMultByNpart->Sumw2();
  hMultByNpartCnt->Sumw2();
  hTRACKp_etHFp  = rpcorrdir.make<TH2D>("TRACKp_etHFp",  "#Psi_{TRACKp} vs.#Psi_{etHFp - TRACKp}",50,-2,2,50,-2,2);
  hTRACKp_HFp    = rpcorrdir.make<TH2D>("TRACKp_HFp",      "#Psi_{TRACKp} vs.#Psi_{HFp - TRACKp}",50,-2,2,50,-2,2);
  hTRACKm_HFm    = rpcorrdir.make<TH2D>("TRACKm_HFm",      "#Psi_{TRACKm} vs.#Psi_{HFm - TRACKm}",50,-2,2,50,-2,2);
  hTRACKm_etHFm  = rpcorrdir.make<TH2D>("TRACKm_etHFm",  "#Psi_{TRACKm} vs. #Psi_{etHFm - TRACKm}",50,-2,2,50,-2,2);
  hTRACKp_TRACKm = rpcorrdir.make<TH2D>("TRACKp_TRACKm","#Psi_{TRACKp} vs. #Psi_{TRACKm - TRACKp}",50,-2,2,50,-2,2);
  hetCaloP_etCaloM = rpcorrdir.make<TH2D>("etCaloP_etCaloM","#Psi_{etCaloP} vs. #Psi_{etCaloP - etCaloM}",50,-2,2,50,-2,2);
  hetHFp_etHFm = rpcorrdir.make<TH2D>("etHFp_etHFm","#Psi_{etHFp} vs. #Psi_{etHFp-etHFm}",50,-2,2,50,-2,2);
  hetHFp_etHFmVsCent = rpcorrdir.make<TH2D>("etHFp_etHFmVsCent","#Psi_{etHFp}-#Psi_{etHFm} vs. cent",nCentBins,centbins,50,-2,2);
  hMid_etHFpVsCent = rpcorrdir.make<TH2D>("Mid_etHFmVsCent","#Psi_{TRACKMid}-#Psi_{etHFp} vs. cent",nCentBins,centbins,50,-2,2);

  hTRACKp_etHFp->SetOption("colz");
  hTRACKp_HFp->SetOption("colz");
  hTRACKm_HFm->SetOption("colz");
  hTRACKm_etHFm->SetOption("colz");
  hTRACKp_TRACKm->SetOption("colz");
  hetCaloP_etCaloM->SetOption("colz");
  hetHFp_etHFm->SetOption("colz");
  hetHFp_etHFmVsCent->SetOption("colz");
  hetHFp_etHFmVsCent->SetXTitle("cent");
  hetHFp_etHFmVsCent->SetYTitle("#Psi_{etHFp}-#Psi_{etHFm}");
  hMid_etHFpVsCent->SetOption("colz");
  hMid_etHFpVsCent->SetXTitle("cent");
  hMid_etHFpVsCent->SetYTitle("#Psi_{TRACKMid}-#Psi_{etHFm}");
  TFileDirectory subdir0 = fs->mkdir("EventPlanes");
  for(int i = 0; i<NumEPNames; i++) {
    for(int j = 0; j<nCentBins; j++) {
      hFlatDiffMean[i][j] = flatchkdir.make<TH1D>(Form("DiffMean_%s_%d_%d",EPNames[i].data(),(int)centbins[j],(int)centbins[j+1]),Form("DiffMean_%s_%d_%d",EPNames[i].data(),(int)centbins[j],(int)centbins[j+1]),100,-2,2);
    }
    hPsi_GenPsi[i] = fs->make<TH2D>(Form("Psi_GenPsi_%s",EPNames[i].data()),Form("Psi_GenPsi_%s",EPNames[i].data()),50,-2,2,50,-2,2); 
    hPsi_GenPsi[i]->SetOption("colz");
    hPsi_GenPsi[i]->SetXTitle(Form("#Psi_{%s}",EPNames[i].data()));
    hPsi_GenPsi[i]->SetYTitle("#Psi_{MC}");
    TFileDirectory subdir = subdir0.mkdir(Form("%s",EPNames[i].data()));
    hFull[i]=subdir.make<TH1D>("psi","psi",200,-4,4);
    hEta[i] = subdir.make<TH1D>("eta","eta",280,-7,7);
    hSub1Sub2[i]=subdir.make<TH2D>("Sub1Sub2","Sub1Sub2",100,-4,4,100,-4,4);
    TFileDirectory subsubdir = subdir.mkdir("CentBins");
    for(int j = 0; j<nCentBins; j++) {
      hFullBin[i][j]=subsubdir.make<TH1D>(Form("psi_%d",j),Form("psi_%d",j),200,-4,4);
      hFullBin[i][j]->Sumw2();
      hSub1Bin[i][j]=subsubdir.make<TH1D>(Form("psiSub1_%d",j),Form("psiSub1_%d",j),200,-4,4);
      hSub1Bin[i][j]->Sumw2();
      hSub2Bin[i][j]=subsubdir.make<TH1D>(Form("psiSub2_%d",j),Form("psiSub2_%d",j),200,-4,4);
      hSub2Bin[i][j]->Sumw2();
    }
    TFileDirectory sub2subdir = subdir.mkdir("GenRes");
    for(int j = 0; j<nCentBins; j++) {
      hGenRes[i][j]=sub2subdir.make<TH1D>(Form("GenRes_%d",j),Form("Gen_%d",j),200,-4,4);
      hGenRes[i][j]->Sumw2();
    }
    TFileDirectory sub3subdir = subdir.mkdir("SubRes");
    for(int j = 0; j<nCentBins; j++) {
      hSubRes[i][j]=sub3subdir.make<TH1D>(Form("SubRes_%d",j),Form("SubRes_%d",j),200,-2,2);
      hSubRes[i][j]->Sumw2();
    }
    TFileDirectory sub4subdir = subdir.mkdir("Mult");
    hMult[i]=sub4subdir.make<TH1D>("Mult","Mult",20,-0.5,19.5);
    hMult[i]->Sumw2();
    hMultCnt[i]=sub4subdir.make<TH1D>("MultCnt","MultCnt",20,-0.5,19.5);
    hMultCnt[i]->Sumw2();
    hMult1[i]=sub4subdir.make<TH1D>("Mult1","Mult1",20,-0.5,19.5);
    hMult1[i]->Sumw2();
    hMult1Cnt[i]=sub4subdir.make<TH1D>("Mult1Cnt","Mult1Cnt",20,-0.5,19.5);
    hMult1Cnt[i]->Sumw2();
    hMult2[i]=sub4subdir.make<TH1D>("Mult2","Mult2",20,-0.5,19.5);
    hMult2[i]->Sumw2();
    hMult2Cnt[i]=sub4subdir.make<TH1D>("Mult2Cnt","Mult2Cnt",20,-0.5,19.5);
    hMult2Cnt[i]->Sumw2();
    TFileDirectory qdir = subdir.mkdir("q");
    for(int j = 0; j<nCentBins; j++) {
      hq[i][j]=qdir.make<TH1D>(Form("q_%d",j),Form("q_%d",j),200,0,5);
      hq[i][j]->Sumw2();
    }
  }
   TFileDirectory v2dir = fs->mkdir("v2");
   TFileDirectory v2Reco = v2dir.mkdir("v2Reco");
   TFileDirectory v2Gen = v2dir.mkdir("v2Gen");
   TFileDirectory PhiEta = v2dir.mkdir("PhiEta");
   TFileDirectory trackv2[NumEPNames] = {
     v2Reco.mkdir(EPNames[0].data()),
     v2Reco.mkdir(EPNames[1].data()),
     v2Reco.mkdir(EPNames[2].data()),
     v2Reco.mkdir(EPNames[3].data()),
     v2Reco.mkdir(EPNames[4].data()),
     v2Reco.mkdir(EPNames[5].data()),
     v2Reco.mkdir(EPNames[6].data()),
     v2Reco.mkdir(EPNames[7].data()),
     v2Reco.mkdir(EPNames[8].data()),
     v2Reco.mkdir(EPNames[9].data()),
     v2Reco.mkdir(EPNames[10].data()),
     v2Reco.mkdir(EPNames[11].data()),
     v2Reco.mkdir(EPNames[12].data()),
     v2Reco.mkdir(EPNames[13].data()),
     v2Reco.mkdir(EPNames[14].data()),
     v2Reco.mkdir(EPNames[15].data()),
     v2Reco.mkdir(EPNames[16].data()),
     v2Reco.mkdir(EPNames[17].data()),
     v2Reco.mkdir(EPNames[18].data()),
     v2Reco.mkdir(EPNames[19].data()),
     v2Reco.mkdir(EPNames[20].data()),
     v2Reco.mkdir(EPNames[21].data()),
     v2Reco.mkdir(EPNames[22].data()),
     v2Reco.mkdir(EPNames[23].data()),
     v2Reco.mkdir(EPNames[24].data()),
     v2Reco.mkdir(EPNames[25].data()),
     v2Reco.mkdir(EPNames[26].data()),
     v2Reco.mkdir(EPNames[27].data()),
     v2Reco.mkdir(EPNames[28].data()),
     v2Reco.mkdir(EPNames[29].data()),
     v2Reco.mkdir(EPNames[30].data()),
     v2Reco.mkdir(EPNames[31].data()),
     v2Reco.mkdir(EPNames[32].data()),
     v2Reco.mkdir(EPNames[33].data()),
     v2Reco.mkdir(EPNames[34].data()),
     v2Reco.mkdir(EPNames[35].data()),
     v2Reco.mkdir(EPNames[36].data()),
     v2Reco.mkdir(EPNames[37].data()),
     v2Reco.mkdir(EPNames[38].data()),
     v2Reco.mkdir(EPNames[39].data()),
     v2Reco.mkdir(EPNames[40].data()),
     v2Reco.mkdir(EPNames[41].data()),
     v2Reco.mkdir(EPNames[42].data()),
     v2Reco.mkdir(EPNames[43].data())
   };
   TFileDirectory calov2[NumEPNames] = {
     v2Reco.mkdir(Form("calo_%s",EPNames[0].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[1].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[2].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[3].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[4].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[5].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[6].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[7].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[8].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[9].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[10].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[11].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[12].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[13].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[14].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[15].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[16].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[17].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[18].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[19].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[20].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[21].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[22].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[23].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[24].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[25].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[26].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[27].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[28].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[29].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[30].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[31].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[32].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[33].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[34].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[35].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[36].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[37].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[38].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[39].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[40].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[41].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[42].data())),
     v2Reco.mkdir(Form("calo_%s",EPNames[43].data()))
   };
   for(Int_t i = 0; i<NumEPNames; i++) {
     v2_Tracks[i]= new v2Generator(trackv2[i],EPNames[i].data(),               0.05);
     v2_Calo[i]  = new v2Generator(calov2[i],Form("calo_%s",EPNames[i].data()),0.05);
   }
   for(Int_t i = 0; i< nPtBins; i++) {
     TString pen = Form("PhiEta_%d",i);
     PhiEta_[i] = PhiEta.make<TH2D>(pen.Data(),pen.Data(),nEtaBins,etabins,180,-TMath::Pi(),TMath::Pi());
     PhiEta_[i]->SetXTitle("#eta");
     PhiEta_[i]->SetYTitle("#phi");
     PhiEta_[i]->SetOption("colz");
   }
}


V2Analyzer::~V2Analyzer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
V2Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace HepMC;
  
  //
  //Get Centrality
  //
  
  if(!centrality_) centrality_ = new CentralityProvider(iSetup);
  centrality_->newEvent(iEvent,iSetup); // make sure you do this first in every event
  //   double c = centrality_->centralityValue();
  int bin = centrality_->getBin();
  //Check if generator events avaiable.  If so, grab info
  const GenEvent *evt;
  Handle<HepMCProduct> mc;
  iEvent.getByLabel("generator",mc);
  double b = 0;
  double npart = 0;
  double Psi = 0;
  double Psi2 = 0;
  double xVert = 0;
  double yVert = 0;
  double zVert = 0;
  double tVert = 0;
  if(mc.isValid()) {
    evt = mc->GetEvent();
    
    const HeavyIon* hi = evt->heavy_ion();
    b = hi->impact_parameter();
    npart = hi->Npart_proj()+hi->Npart_targ();
    Psi = hi->event_plane_angle();
    Psi2 = bounds2(Psi);
    HepMC::GenVertex * sigvert = evt->signal_process_vertex();
    if(sigvert) {
      HepMC::FourVector point = sigvert->position();
      xVert=point.x();
      yVert=point.y();
      zVert=point.z();
      tVert=point.t();
    }
  }
  
  double centval = 2.5*bin+1.25;
  hNpartBin->Fill( centval,centrality_->NpartMean() );
  hNpartBinCnt->Fill(centval);
  bin = hCentBinned->FindBin(centval)-1;
  if(bin<0) return;
  hcent->Fill(centval);
  //
  //Get Vertex
  //
  edm::Handle<reco::VertexCollection> vertexCollection3;
  iEvent.getByLabel("hiSelectedVertex",vertexCollection3);
  const reco::VertexCollection * vertices3 = vertexCollection3.product();
  vs_sell = vertices3->size();
  if(vs_sell>0) {
    vzr_sell = vertices3->begin()->z();
    vzErr_sell = vertices3->begin()->zError();
  } else
    vzr_sell = -999.9;
  //
  //Get Event Planes
  //
  Handle<reco::EvtPlaneCollection> evtPlanes;
  Handle<reco::EvtPlaneCollection> evtPlanesNoFlat;
  iEvent.getByLabel("hiEvtPlaneFlat","recoLevel",evtPlanes);
  iEvent.getByLabel("hiEvtPlane","recoLevel",evtPlanesNoFlat);
  double fullNoFlat[NumEPNames];
  if(evtPlanesNoFlat.isValid()){
    for (EvtPlaneCollection::const_iterator rp = evtPlanesNoFlat->begin();rp !=evtPlanesNoFlat->end(); rp++) {
      size_t pos;
      if(rp->angle() > -5) {
	pos = rp->label().find("_sub");      
	string baseName = rp->label();
	if(pos != string::npos) baseName = rp->label().substr(0,pos);
	for(int i = 0; i< NumEPNames; i++) {
	  if(EPNames[i].compare(baseName)==0) {
	    if(EPNames[i].compare(rp->label())==0) {
	      fullNoFlat[i]=rp->angle();
	    }
	  }
	}    
      }
    }
  }
  if(!evtPlanes.isValid()){
    cout << "Error! Can't get hiEvtPlane product!" << endl;
    return ;
  }
  double full[NumEPNames];
  double sub1[NumEPNames];
  double sub2[NumEPNames];
  double mult[NumEPNames];
  double mult1[NumEPNames];
  double mult2[NumEPNames];
  double sumSin[NumEPNames];
  double sumCos[NumEPNames];
  double Q[NumEPNames];
  
  for(int i = 0; i<NumEPNames;i++) {
    full[i] = -10;
    sub1[i]=-10;
    sub2[i]=-10;
    mult[i]=0;
    mult1[i]=0;
    mult2[i]=0;
  }
  for (EvtPlaneCollection::const_iterator rp = evtPlanes->begin();rp !=evtPlanes->end(); rp++) {
    size_t pos;
    if(rp->angle() > -5) {
      pos = rp->label().find("_sub");      
      string baseName = rp->label();
      if(pos != string::npos) baseName = rp->label().substr(0,pos);
      for(int i = 0; i< NumEPNames; i++) {
	if(EPNames[i].compare(baseName)==0) {
	  double multScale = 1.;
	  if(baseName.find("et")!=string::npos) multScale=1./0.9;  //convert calo et to mult
	  if(EPNames[i].compare(rp->label())==0) {
	    full[i]=rp->angle();
	    mult[i]=rp->mult()*multScale;
	    sumSin[i]=rp->sumSin();
	    sumCos[i]=rp->sumCos();
	    Q[i]=rp->Q();
	    if(mc.isValid()) hPsi_GenPsi[i]->Fill(full[i],Psi2);
	  } else if (rp->label().find("_sub1") != string::npos) {
	    sub1[i]=rp->angle();
	    mult1[i]=rp->mult()*multScale;
	  } else if (rp->label().find("_sub2") != string::npos) {
	    sub2[i]=rp->angle();
	    mult2[i]=rp->mult()*multScale;
	  }
	}
      }
    }    
  }
  for(int i = 0; i< 8; i++) {
    hPsiFlatVsPsiRaw[i]->Fill(fullNoFlat[i],full[i]);
  }
  if (centval>=50 && centval<60) {
    hTRACKp_etHFp->Fill(bounds2(full[etHFp]-full[EvtPTracksPosEtaGap]),full[EvtPTracksPosEtaGap]);
    hTRACKp_HFp->Fill(bounds2(full[HFp]-full[EvtPTracksPosEtaGap]),full[EvtPTracksPosEtaGap]);
    hTRACKm_HFm->Fill(bounds2(full[HFm]-full[EvtPTracksNegEtaGap]),full[EvtPTracksNegEtaGap]);
    hTRACKm_etHFm->Fill(bounds2(full[etHFm]-full[EvtPTracksNegEtaGap]),full[EvtPTracksNegEtaGap]);
    hTRACKp_TRACKm->Fill(bounds2(full[EvtPTracksNegEtaGap]-full[EvtPTracksPosEtaGap]),full[EvtPTracksPosEtaGap]);
  }
  if(full[etHFp]>-5&&full[etHFm]>-5) hetHFp_etHFmVsCent->Fill(centval,bounds2(full[etHFp]-full[etHFm]));
  if(full[EvtPlaneFromTracksMidEta]>-5&&full[etHFm]>-5) 
    hMid_etHFpVsCent->Fill(centval,bounds2(full[EvtPlaneFromTracksMidEta]-full[etHFm]));
  if(centval>5&&centval<40) {
    hetCaloP_etCaloM->Fill(bounds2(full[etCaloP]-full[etCaloM]),bounds2(full[etCaloP]));
    hetHFp_etHFm->Fill(bounds2(full[etHFp]-full[etHFm]),bounds2(full[etHFp]));
  }
  for(int i = 0; i< NumEPNames; i++) {
    if(full[i]>-5) { 
      Double_t diff = full[i]-fullNoFlat[i];
      hFlatDiffMean[i][bin]->Fill(diff);
  
      hFull[i]->Fill(full[i]);
      if(sub1[i]>-5 && sub2[i]>-5) hSub1Sub2[i]->Fill(sub1[i]-Psi2,sub2[i]-Psi2);
      hFullBin[i][bin]->Fill(full[i]);
      hMult[i]->Fill(bin,mult[i]);
      hMultCnt[i]->Fill(bin);
      if(mult[i]>0) hq[i][bin]->Fill(Q[i]/mult[i]);
      if(EPNames[i].find("1")!=string::npos) {
	hGenRes[i][bin]->Fill(cos(full[i]-Psi));
      } else {
	hGenRes[i][bin]->Fill(cos(2.*(full[i]-Psi2)));
      }
      if(sub1[i]>-5 && sub2[i]>-5) {
	hSubRes[i][bin]->Fill(cos(2.*(sub1[i] - sub2[i])));
	hMult1[i]->Fill(bin,mult1[i]);
	hMult2[i]->Fill(bin,mult2[i]);
	hMult1Cnt[i]->Fill(bin);
	hMult2Cnt[i]->Fill(bin);
	hSub1Bin[i][bin]->Fill(sub1[i]);
	hSub2Bin[i][bin]->Fill(sub2[i]);
      }
      v2_Tracks[i]->AddToResCor(full[i],full[RCMate1[i]],full[RCMate2[i]],centval);
      v2_Calo[i]->AddToResCor(full[i],full[RCMate1[i]],full[RCMate2[i]],centval);
      if(mc.isValid()) {
	v2_Tracks[i]->AddToGenRes(full[i],Psi2,centval);
	v2_Calo[i]->AddToGenRes(full[i],Psi2,centval);
      }
    }
  }
  //Tracking part
  double track_eta=-10;
  double track_phi=-10;
  double track_pt=-10;
  //double track_charge;
  
  // for(int i = 0; i<8; i++) {
    //v2_Tracks[i]->ResetAutocorrelation();
    //v2_etCaloHF[i]->ResetAutocorrelation();
  // }
#ifdef TRACKCOLLECTION  
  Handle<reco::TrackCollection> tracks;
  //  iEvent.getByLabel("hiSelectedTracks", tracks);
  iEvent.getByLabel("hiGoodMergedTracks", tracks);
  if(tracks.isValid()){
    //for(reco::TrackCollection::const_iterator k = tracks->begin(); k!= tracks->end(); k++) {
      //      double w = 1;
      //for(int i = 0; i< NumEPNames; i++) v2_Tracks[i]->SetAutocorrelation(k->phi(), k->eta(), w);
    //}
    for(reco::TrackCollection::const_iterator j = tracks->begin(); j != tracks->end(); j++){
#endif
#ifdef RECOCHARGEDCANDIDATECOLLECTION
   edm::Handle<reco::RecoChargedCandidateCollection> trackCollection;
   iEvent.getByLabel("allMergedPtSplit12Tracks",trackCollection);
   
   if(trackCollection.isValid()){
     const reco::RecoChargedCandidateCollection * tracks = trackCollection.product();
     //for(reco::RecoChargedCandidateCollection::const_iterator k = tracks->begin(); k!= tracks->end(); k++) {
     //  double w = 1;
     //  for(int i = 0; i< NumEPNames; i++) v2_Tracks[i]->SetAutocorrelation(k->phi(), k->eta(), w);
     //}
     for(reco::RecoChargedCandidateCollection::const_iterator j = tracks->begin(); j != tracks->end(); j++){
#endif     
       edm::Handle<reco::VertexCollection> vertex;
       iEvent.getByLabel("hiSelectedVertex", vertex);

// find the vertex point and error

       math::XYZPoint vtxPoint(0.0,0.0,0.0);
       double vzErr =0.0, vxErr=0.0, vyErr=0.0;
       if(vertex->size()>0) {
	 vtxPoint=vertex->begin()->position();
	 vzErr=vertex->begin()->zError();
	 vxErr=vertex->begin()->xError();
	 vyErr=vertex->begin()->yError();
       }
       bool accepted = true;
       bool isPixel = false;
       // determine if the track is a pixel track
       if ( j->numberOfValidHits() < 7 ) isPixel = true;
       
       // determine the vertex significance 
       double d0=0.0, dz=0.0, d0sigma=0.0, dzsigma=0.0;
       d0 = -1.*j->dxy(vtxPoint);
       dz = j->dz(vtxPoint);
       d0sigma = sqrt(j->d0Error()*j->d0Error()+vxErr*vyErr);
       dzsigma = sqrt(j->dzError()*j->dzError()+vzErr*vzErr);
       
       // cuts for pixel tracks
       if( isPixel )
	 {
	   // dz significance cut 
	   if ( fabs(dz/dzsigma) > dzerr_ ) accepted = false;
	   
	   // chi2/ndof cut 
	   if ( j->normalizedChi2() > chi2_ ) accepted = false;
	 }
       
       // cuts for full tracks
       if ( ! isPixel)
	 {
	   // dz and d0 significance cuts 
	   if ( fabs(dz/dzsigma) > 3 ) accepted = false;
	   if ( fabs(d0/d0sigma) > 3 ) accepted = false;
	   
	   // pt resolution cut
	   if ( j->ptError()/j->pt() > 0.05 ) accepted = false;
	   
	   // number of valid hits cut
	   if ( j->numberOfValidHits() < 12 ) accepted = false;
	 }
       
       if( accepted ){     
	 
	 track_eta = j->eta();
	 track_phi = j->phi();
	 track_pt = j->pt();
	 int ptbin = hpt[0]->GetXaxis()->FindBin(track_pt)-1;
	 if(ptbin>=0&&ptbin<nPtBins) PhiEta_[ptbin]->Fill(track_eta,track_phi);
	 cout<<"valid track: "<<track_eta<<" "<<track_phi<<" "<<track_pt<<endl;
	 //track_charge = j->charge();
	 for(int i = 0; i< NumEPNames; i++) {
	   TString ename = EPNames[i].data();
	   if(ename.Contains("Track") && InEtaRange(i,track_eta)) hEta[i]->Fill(track_eta);
	   if(track_eta>1.8&&track_phi>-0.1&&track_phi<1.0) continue;
	   v2_Tracks[i]->AddParticle(track_phi,full[i],centval,track_eta,track_pt);
	   if(mc.isValid()) {
	     v2_Tracks[i]->AddGenParticle(track_phi,Psi2,centval,track_eta,track_pt);
	   }
	 }
	 Int_t ietabin = heta->FindBin(track_eta)-1;
	 if(ietabin>=0) {
	   hpt[ietabin]->Fill(track_pt,centval,track_pt);
	   hptCnt[ietabin]->Fill(track_pt,centval);
	 }
       }
     }
   } else {
     cout<<"Failed to find trackCollection"<<endl;
   }
   
  
  Handle<CaloTowerCollection> calotower;
  iEvent.getByLabel("towerMaker",calotower);
  if(calotower.isValid()){
    for (CaloTowerCollection::const_iterator j = calotower->begin();j !=calotower->end(); j++) {   
      double w = j->emEt()+j->hadEt();
      for(int i = 0; i<NumEPNames; i++) {
	TString ename = EPNames[i].data();
	if(!ename.Contains("Track") && InEtaRange(i,track_eta)) hEta[i]->Fill(j->eta());
	v2_Calo[i]->AddParticle(j->phi(),full[i],centval,j->eta(),w);
	if(mc.isValid()) {
	  v2_Calo[i]->AddGenParticle(j->phi(),Psi2,centval,j->eta(),w);
	}
      }
      Int_t ietabin = heta->FindBin(j->eta())-1;
      if(ietabin>=0) {
       	het[ietabin]->Fill(w,centval,w);
       	hetCnt[ietabin]->Fill(w,centval);
	hEmHad_EmEt[ietabin]->Fill(w,j->emEt());
      }
    }
    
  }
} 

  // ------------ method called once each job just before starting event loop  ------------
  void 
    V2Analyzer::beginJob()
  {
  }
  
  // ------------ method called once each job just after ending the event loop  ------------
  void 
    V2Analyzer::endJob() {
  }

  //define this as a plug-in
  DEFINE_FWK_MODULE(V2Analyzer);
  
