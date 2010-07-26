// -*- C++ -*-
//
// Package:    V2Analyzer
// Class:      V2Analyzer
// 
/**\class V2Analyzer V2Analyzer.cc HiEvtPlaneFlatten/V2Analyzer/src/V2Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Stephen Sanders
//         Created:  Wed Jun 23 12:27:13 EDT 2010
// $Id: V2Analyzer.cc,v 1.1 2010/07/19 22:11:12 ssanders Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
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

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"

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
  
  class v2Generator {
  public:
    explicit v2Generator(TFileDirectory  dir, string label, double gap,  int nCentBins, double * centBins, int nEtaBins, double * etaBins, int nPtBins, double * ptBins){
      label_ = label;
      gap_ = fabs(gap);
      if(gap_ == 0) {
	subAuto = false;
      } else {
	subAuto = true;
	static const int maxAutoEtaBins = 4000;
	int numAutoEtaBins = (etaBins[nEtaBins]-etaBins[0])/gap_;
	if(numAutoEtaBins > maxAutoEtaBins) {
	  gap_ = (etaBins[nEtaBins]-etaBins[0])/maxAutoEtaBins;
	  numAutoEtaBins = maxAutoEtaBins;
	}
	autoSin = dir.make<TH1D>("autoSin","autoSin",4*numAutoEtaBins,etaBins[0],etaBins[nEtaBins]);
	autoSin->Sumw2();
	autoCos = dir.make<TH1D>("autoCos","autoCos",4*numAutoEtaBins,etaBins[0],etaBins[nEtaBins]);
	autoCos->Sumw2();
 	autoCnt = dir.make<TH1D>("autoCnt","autoCnt",4*numAutoEtaBins,etaBins[0],etaBins[nEtaBins]);
	autoCnt->Sumw2();
	autoPsi = dir.make<TH1D>("autoPsi","Psi distribution after autocorrelation correction",200,-4,4);
	autoPsi->Sumw2();
     }
      nEtaBins_ = nEtaBins; 
      eta_ = dir.make<TH1F>(Form("eta_%s",label.data()),Form("eta_%s",label.data()),nEtaBins,etaBins);
      for(int i = 0; i<nEtaBins; i++) {
	cos_[i] = dir.make<TH2D>(Form("cos_%s_%d",label.data(),i),Form("cos_%s_%d",label.data(),i),nCentBins,centBins,nPtBins,ptBins);
	cos_[i]->SetXTitle("cent");
	cos_[i]->SetYTitle("p_{T}");
	cos_[i]->SetOption("colz");
	cos_[i]->Sumw2();
	cnt_[i] = dir.make<TH2D>(Form("cnt_%s_%d",label.data(),i),Form("cnt_%s_%d",label.data(),i),nCentBins,centBins,nPtBins,ptBins);
	cnt_[i]->SetXTitle("cent");
	cnt_[i]->SetYTitle("p_{T}");
	cnt_[i]->SetOption("colz");
	cnt_[i]->Sumw2();
      } 
    }
    ~v2Generator();
   
    void AddParticle(double phi, double Psi, double cent, double eta, double pt) {
      int ietabin = eta_->FindBin(eta)-1;
      if(ietabin<0 || ietabin>=nEtaBins_) return;
      cos_[ietabin]->Fill(cent,pt,TMath::Cos(2*(phi-Psi)));
      cnt_[ietabin]->Fill(cent,pt);
    }
    void SetAutocorrelation(double phi, double eta, double w) {
      if(!subAuto) return;
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
    string label_;
    double gap_;
    TH2D * cos_[50];
    TH2D * cnt_[50];
    int nEtaBins_; 
    TH1F * eta_;
    TH1D * autoSin;
    TH1D * autoCos;
    TH1D * autoCnt;
    TH1D * autoPsi; 
    bool subAuto;
  };
  
  edm::Service<TFileService> fs;
  const CentralityBins * cbins_;
  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TH1D * hcent;
  TH1D * hMultByNpart;
  TH1D * hMultByNpartCnt;
  TH1D * hFull[NumEPNames];
  TH1D * hFullBin[NumEPNames][20];
  TH1D * hSub1Bin[NumEPNames][20];
  TH1D * hSub2Bin[NumEPNames][20];
  TH2D * hSub1Sub2[NumEPNames];
  TH1D * hGenRes[NumEPNames][20];
  TH1D * hSubRes[NumEPNames][20];
  TH1D * hMult1[NumEPNames];
  TH1D * hMult2[NumEPNames];
  TH1D * hMult[NumEPNames];
  TH1D * hMult1Cnt[NumEPNames];
  TH1D * hMult2Cnt[NumEPNames];
  TH1D * hMultCnt[NumEPNames];
  TH1D * hq[NumEPNames][20];
  TH2D * hEPCorrelation[20];
  TH2D * hEPCorrelationCnt[20];
  TH2D * hv2RecoCos;
  TH2D * hv2RecoCnt;
  TH2D * hv2GenCos;
  TH2D * hv2GenCnt;
  TH2D * hv1RecoCos;
  TH2D * hv1RecoCnt;
  TH2D * hv1GenCos;
  TH2D * hv1GenCnt;
  TH1D * hNpartBin;
  TH1D * hNpartBinCnt;
  TH2D * hpt;
  TH2D * hptCnt;
  v2Generator * v2_Tracks[8];
  v2Generator * v2_etCaloHF;
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
   //now do what ever initialization is needed

  double centbins[]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};
  double ptbins[]={0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.6,2.0,2.5,3.0,4.0,6.0,8.0,12.0};
  double etabinsTracks[]={-2,-1,0,1,2};
  cbins_ = 0;
  hcent = fs->make<TH1D>("cent","cent",200,-10,110);
  hNpartBin = fs->make<TH1D>("NpartBin","NpartBin",41,-0.5,40.5);
  hNpartBinCnt = fs->make<TH1D>("NpartBinCnt","NpartBinCnt",41,-0.5,40.5);
  hpt = fs->make<TH2D>("pt","pt",15,ptbins,20,centbins);
  hptCnt = fs->make<TH2D>("ptCnt","ptCnt",15,ptbins,20,centbins);
  hpt->Sumw2();
  hptCnt->Sumw2();
  hMultByNpart = fs->make<TH1D>("MultByNpart","2*Mult/N_{part}",200,-10,110);
  hMultByNpartCnt = fs->make<TH1D>("MultByNpartCnt","MultByNpartCnt",200,-10,110);
  TFileDirectory subdir0 = fs->mkdir("EventPlanes");
  for(int i = 0; i<NumEPNames; i++) {
    TFileDirectory subdir = subdir0.mkdir(Form("%s",EPNames[i].data()));
    hFull[i]=subdir.make<TH1D>("psi","psi",200,-4,4);
    hSub1Sub2[i]=subdir.make<TH2D>("Sub1Sub2","Sub1Sub2",100,-4,4,100,-4,4);
    TFileDirectory subsubdir = subdir.mkdir("CentBins");
    for(int j = 0; j<20; j++) {
      hFullBin[i][j]=subsubdir.make<TH1D>(Form("psi_%d",j),Form("psi_%d",j),200,-4,4);
      hFullBin[i][j]->Sumw2();
      hSub1Bin[i][j]=subsubdir.make<TH1D>(Form("psiSub1_%d",j),Form("psiSub1_%d",j),200,-4,4);
      hSub1Bin[i][j]->Sumw2();
      hSub2Bin[i][j]=subsubdir.make<TH1D>(Form("psiSub2_%d",j),Form("psiSub2_%d",j),200,-4,4);
      hSub2Bin[i][j]->Sumw2();
    }
    TFileDirectory sub2subdir = subdir.mkdir("GenRes");
    for(int j = 0; j<20; j++) {
      hGenRes[i][j]=sub2subdir.make<TH1D>(Form("GenRes_%d",j),Form("Gen_%d",j),200,-4,4);
      hGenRes[i][j]->Sumw2();
    }
    TFileDirectory sub3subdir = subdir.mkdir("SubRes");
    for(int j = 0; j<20; j++) {
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
    for(int j = 0; j<20; j++) {
      hq[i][j]=qdir.make<TH1D>(Form("q_%d",j),Form("q_%d",j),200,0,5);
      hq[i][j]->Sumw2();
    }
  }
  TFileDirectory subdir = fs->mkdir("EPCorrelations");
  for(int j = 0; j<20; j++) {
    hEPCorrelation[j] = subdir.make<TH2D>(Form("EPCorrelation_%d",j),Form("EPCorrelation_%d",j),NumEPNames,-0.5,NumEPNames-0.5,NumEPNames,-0.5,NumEPNames-0.5);
    hEPCorrelationCnt[j] = subdir.make<TH2D>(Form("EPCorrelationCnt_%d",j),Form("EPCorrelationCnt_%d",j),NumEPNames,-0.5,NumEPNames-0.5,NumEPNames,-0.5,NumEPNames-0.5);
    hEPCorrelation[j]->SetStats(kFALSE);
    hEPCorrelation[j]->SetOption("colz");
    hEPCorrelation[j]->Sumw2();
    hEPCorrelationCnt[j]->Sumw2();
  }
  TFileDirectory v2dir = fs->mkdir("v2");
  TFileDirectory v2Reco = v2dir.mkdir("v2Reco");
  TFileDirectory v2Gen = v2dir.mkdir("v2Gen");

  hv2RecoCos = v2Reco.make<TH2D>("v2RecoCos","v2RecoCos",20,centbins,4,etabinsTracks);
  hv2RecoCos->Sumw2();
  hv2RecoCnt = v2Reco.make<TH2D>("v2RecoCnt","v2RecoCnt",20,centbins,4,etabinsTracks);
  hv2RecoCnt->Sumw2();
  hv2GenCos = v2Reco.make<TH2D>("v2GenCos","v2GenCos",20,centbins,4,etabinsTracks);
  hv2GenCos->Sumw2();
  hv2GenCnt = v2Reco.make<TH2D>("v2GenCnt","v2GenCnt",20,centbins,4,etabinsTracks);
  hv2GenCnt->Sumw2();

  TFileDirectory trackv2_0 = v2Reco.mkdir(EPNames[0].data()); 
  TFileDirectory trackv2_1 = v2Reco.mkdir(EPNames[1].data()); 
  TFileDirectory trackv2_2 = v2Reco.mkdir(EPNames[2].data()); 
  TFileDirectory trackv2_3 = v2Reco.mkdir(EPNames[3].data()); 
  TFileDirectory trackv2_4 = v2Reco.mkdir(EPNames[4].data()); 
  TFileDirectory trackv2_5 = v2Reco.mkdir(EPNames[5].data()); 
  TFileDirectory trackv2_6 = v2Reco.mkdir(EPNames[6].data()); 
  TFileDirectory trackv2_7 = v2Reco.mkdir(EPNames[7].data()); 
  v2_Tracks[0] = new v2Generator(trackv2_0,EPNames[0].data(),0.05,20,centbins,4,etabinsTracks,15,ptbins);
  v2_Tracks[1] = new v2Generator(trackv2_1,EPNames[1].data(),0.05,20,centbins,4,etabinsTracks,15,ptbins);
  v2_Tracks[2] = new v2Generator(trackv2_2,EPNames[2].data(),0.05,20,centbins,4,etabinsTracks,15,ptbins);
  v2_Tracks[3] = new v2Generator(trackv2_3,EPNames[3].data(),0.05,20,centbins,4,etabinsTracks,15,ptbins);
  v2_Tracks[4] = new v2Generator(trackv2_4,EPNames[4].data(),0.05,20,centbins,4,etabinsTracks,15,ptbins);
  v2_Tracks[5] = new v2Generator(trackv2_5,EPNames[5].data(),0.05,20,centbins,4,etabinsTracks,15,ptbins);
  v2_Tracks[6] = new v2Generator(trackv2_6,EPNames[6].data(),0.05,20,centbins,4,etabinsTracks,15,ptbins);
  v2_Tracks[7] = new v2Generator(trackv2_7,EPNames[7].data(),0.05,20,centbins,4,etabinsTracks,15,ptbins);

  TFileDirectory caloHFv2 = v2Reco.mkdir(EPNames[etCaloHF].data());
  v2_etCaloHF = new v2Generator(caloHFv2,EPNames[etCaloHF].data(),0.1,20,centbins,4,etabinsTracks,15,ptbins);


  TFileDirectory v1dir = fs->mkdir("v1");
  TFileDirectory v1Reco = v1dir.mkdir("v1Reco");
  TFileDirectory v1Gen = v1dir.mkdir("v1Gen");

  hv1RecoCos = v1Reco.make<TH2D>("v1RecoCos","v1RecoCos",20,centbins,4,etabinsTracks);
  hv1RecoCos->Sumw2();
  hv1RecoCnt = v1Reco.make<TH2D>("v1RecoCnt","v1RecoCnt",20,centbins,4,etabinsTracks);
  hv1RecoCnt->Sumw2();
  hv1GenCos = v1Reco.make<TH2D>("v1GenCos","v1GenCos",20,centbins,4,etabinsTracks);
  hv1GenCos->Sumw2();
  hv1GenCnt = v1Reco.make<TH2D>("v1GenCnt","v1GenCnt",20,centbins,4,etabinsTracks);
  hv1GenCnt->Sumw2();
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
   //
  //Get Centrality
  //
   if(!cbins_) cbins_ = getCentralityBinsFromDB(iSetup);
   edm::Handle<reco::Centrality> cent;
   iEvent.getByLabel(edm::InputTag("hiCentrality"),cent);
   if(!cent.isValid()){
     cout << "Error! Can't get hiCentrality product!" << endl;
     return ;
   }  
   double  hf = cent->EtHFhitSum();
   int bin = cbins_->getBin(hf);
   hNpartBin->Fill( bin,cbins_->NpartMeanOfBin(bin) );
   hNpartBinCnt->Fill(bin);
   double centval = 5.*bin+0.25;
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
  iEvent.getByLabel("hiEvtPlaneFlat","recoLevel",evtPlanes);
  //iEvent.getByLabel("hiEvtPlane","recoLevel",evtPlanes);

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

  for(int i = 0; i< NumEPNames; i++) {
    if(full[i]>-5) { 
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
    }
  }
  for(int i = 0; i< NumEPNames; i++) {
    if(EPNames[i].find("1")==string::npos && full[i]>-5){
      for(int j = 0; j<NumEPNames; j++) {
	if(EPNames[j].find("1") == string::npos && full[j]>-5){
	  hEPCorrelation[bin]->Fill(i,j,cos(2.*(full[i]-full[j])));
	  hEPCorrelationCnt[bin]->Fill(i,j);
	}
      }
    }
  }
  //Tracking part
  
  double track_eta;
  double track_phi;
  double track_pt;
  //double track_charge;
  
  for(int i = 0; i<8; i++) v2_Tracks[i]->ResetAutocorrelation();
  v2_etCaloHF->ResetAutocorrelation();
  
  //Handle<reco::TrackCollection> tracks;
  //iEvent.getByLabel("hiSelectedTracks", tracks);

  edm::Handle<reco::RecoChargedCandidateCollection> trackCollection;
  iEvent.getByLabel("allMergedPtSplit12Tracks",trackCollection);

  //  if(tracks.isValid()){
  if(trackCollection.isValid()){
    const reco::RecoChargedCandidateCollection * tracks = trackCollection.product();
    //    for(reco::TrackCollection::const_iterator k = tracks->begin(); k!= tracks->end(); k++) {
    for(reco::RecoChargedCandidateCollection::const_iterator k = tracks->begin(); k!= tracks->end(); k++) {
      double w = 1;
      for(int i = 0; i< 8; i++) v2_Tracks[i]->SetAutocorrelation(k->phi(), k->eta(), w);
    }
    //    for(reco::TrackCollection::const_iterator j = tracks->begin(); j != tracks->end(); j++){
    for(reco::RecoChargedCandidateCollection::const_iterator j = tracks->begin(); j != tracks->end(); j++){
      
      track_eta = j->eta();
      track_phi = j->phi();
      track_pt = j->pt();
      //track_charge = j->charge();
      double psiReco = -10;
      int trackbin = (int) 2*(track_eta+2);
      if(trackbin>=0 && trackbin < 8) {
	hv1GenCos->Fill(5*bin+2.5,track_eta,cos(track_phi - Psi2));
	hv1GenCnt->Fill(5*bin+2.5,track_eta);
	hv2GenCos->Fill(5*bin+2.5,track_eta,cos(2.0*(track_phi - Psi2)));
	hv2GenCnt->Fill(5*bin+2.5,track_eta);
	if(trackbin>=0 && trackbin<4){
	  psiReco = full[ EvtPTracksPosEtaGap ];
	  hMultByNpart->Fill(centval, mult[EvtPlaneFromTracksEta]/cbins_->NpartMeanOfBin(bin));
	  hMultByNpartCnt->Fill(centval); 
	} else if (trackbin >=4 && trackbin < 8) {
	  psiReco = full[ EvtPTracksNegEtaGap ];
	}
	//	if(psiReco > -5) {
	//	  hv2RecoCos->Fill(5*bin+2.5,track_eta,cos(2.0*(track_phi - psiReco)));
	//	  hv2RecoCnt->Fill(5*bin+2.5,track_eta);
	//	  hv1RecoCos->Fill(5*bin+2.5,track_eta,cos(track_phi - psiReco));
	//	  hv1RecoCnt->Fill(5*bin+2.5,track_eta);
	//}
      }
      for(int i = 0; i< 8; i++) {
	psiReco = full[i];
	//psiReco = v2_Tracks[i]->GetAutoCorrectedPsi(track_eta, full[i], sumSin[i], sumCos[i]);
	v2_Tracks[i]->AddParticle(track_phi,psiReco,5*bin+2.5,track_eta,track_pt);
      }
      if(fabs(track_eta)<1.) {
	hpt->Fill(track_pt,5*bin+2.5,track_pt);
	hptCnt->Fill(track_pt,5*bin+2.5);
      }
      psiReco = v2_etCaloHF->GetAutoCorrectedPsi(track_eta, full[etCaloHF], sumSin[etCaloHF], sumCos[etCaloHF]);
      v2_etCaloHF->AddParticle(track_phi,psiReco,5*bin+2.5,track_eta,track_pt);
    }
  }
  
  Handle<CaloTowerCollection> calotower;
  iEvent.getByLabel("towerMaker",calotower);
  if(calotower.isValid()){
    for (CaloTowerCollection::const_iterator j = calotower->begin();j !=calotower->end(); j++) {   
      double w = j->emEt()+j->hadEt();
      v2_etCaloHF->SetAutocorrelation(j->phi(), j->eta(), w);
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
