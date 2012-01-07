#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TVector3.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Utilities/General/interface/FileInPath.h"
#include "HiggsAnalysis/HiggsToWW2Leptons/interface/ElectronIDMVA.h"

#include "VHTauTau/TreeMaker/plugins/ElectronBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

// Constructor
ElectronBlock::ElectronBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _trkInputTag(iConfig.getParameter<edm::InputTag>("trackSrc")),
  _dcsInputTag(iConfig.getParameter<edm::InputTag>("dcsSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  _electronInputTag(iConfig.getParameter<edm::InputTag>("electronSrc")),
  _pfElectronInputTag(iConfig.getParameter<edm::InputTag>("pfElectronSrc")),
  _ecalEBInputTag(iConfig.getParameter<edm::InputTag>("ecalEBInputTag")),
  _ecalEEInputTag(iConfig.getParameter<edm::InputTag>("ecalEEInputTag"))
{
  std::string method = iConfig.getParameter<std::string>("methodName");
  edm::FileInPath Subdet0Pt10To20Weights = iConfig.getParameter<edm::FileInPath>("Subdet0LowPtWeights");
  edm::FileInPath Subdet1Pt10To20Weights = iConfig.getParameter<edm::FileInPath>("Subdet1LowPtWeights");
  edm::FileInPath Subdet2Pt10To20Weights = iConfig.getParameter<edm::FileInPath>("Subdet2LowPtWeights");
  edm::FileInPath Subdet0HighPtWeights = iConfig.getParameter<edm::FileInPath>("Subdet0HighPtWeights");
  edm::FileInPath Subdet1HighPtWeights = iConfig.getParameter<edm::FileInPath>("Subdet1HighPtWeights");
  edm::FileInPath Subdet2HighPtWeights = iConfig.getParameter<edm::FileInPath>("Subdet2HighPtWeights");
  ElectronIDMVA::MVAType mvaType = static_cast<ElectronIDMVA::MVAType>(iConfig.getParameter<unsigned int>("mvaType")); 

  // Electron ID MVA against fake QCD
  fMVA = new ElectronIDMVA();
  fMVA->Initialize(method,
		   Subdet0Pt10To20Weights.fullPath(),
		   Subdet1Pt10To20Weights.fullPath(),
		   Subdet2Pt10To20Weights.fullPath(),
		   Subdet0HighPtWeights.fullPath(),
		   Subdet1HighPtWeights.fullPath(),
		   Subdet2HighPtWeights.fullPath(),
                   mvaType);
}
ElectronBlock::~ElectronBlock() { delete fMVA; }
void ElectronBlock::beginJob() 
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  cloneElectron = new TClonesArray("vhtm::Electron");
  tree->Branch("Electron", &cloneElectron, 32000, 2);
  tree->Branch("nElectron", &fnElectron,  "fnElectron/I");
}
void ElectronBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneElectron->Clear();
  fnElectron = 0;

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(_trkInputTag, tracks);

  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByLabel(_dcsInputTag, dcsHandle);

  double evt_bField = 3.8;
  // need the magnetic field
  //
  // if isRealData then derive bfield using the
  // magnet current from DcsStatus
  // otherwise take it from the IdealMagneticFieldRecord
  if (iEvent.isRealData()) {
    if (dcsHandle.isValid()) {
      edm::LogInfo("ElectronBlock") << "Successfully obtained " << _dcsInputTag;
      // scale factor = 3.801/18166.0 which are
      // average values taken over a stable two-week period
      double currentToBFieldScaleFactor = 2.09237036221512717e-04;
      double current = (*dcsHandle)[0].magnetCurrent();
      evt_bField = current*currentToBFieldScaleFactor;
    } 
    else {
      edm::LogError("ElectronBlock") << "Error! Can't get the product " << _dcsInputTag;
    }
  } 
  else {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    if (magneticField.isValid()) {
      edm::LogInfo("ElectronBlock") << "Successfully obtained IdealMagneticFieldRecord!";
      evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
    } 
    else {
      edm::LogError("ElectronBlock") << "Error! Can't get IdealMagneticFieldRecord";
    }
  }

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(_electronInputTag, electrons);

  edm::Handle<reco::PFCandidateCollection> pfElectrons;
  iEvent.getByLabel(_pfElectronInputTag, pfElectrons);

  if (electrons.isValid()) {
    edm::LogInfo("ElectronBlock") << "Total # PAT Electrons: " << electrons->size();
    for (std::vector<pat::Electron>::const_iterator it = electrons->begin(); 
                                                   it != electrons->end(); ++it) {
      if (fnElectron == kMaxElectron) {
	edm::LogInfo("ElectronBlock") << "Too many PAT Electrons, fnElectron = " << fnElectron; 
	break;
      }
      // if electron is not ECAL driven, continue
      if (!it->ecalDrivenSeed()) continue;

      // Conversion
      ConversionFinder convFinder;
      double dist = -9999.;
      double dcot = -9999.;
      if (tracks.isValid()) {
	edm::LogInfo("ElectronBlock") << "Successfully obtained " << _trkInputTag;

        ConversionInfo convInfo = convFinder.getConversionInfo(*it, tracks, evt_bField);
        dist = convInfo.dist();
        dcot = convInfo.dcot();
      } 
      else {
	edm::LogError("ElectronBlock") << "Error! Can't get the product " << _trkInputTag;
      }

      // Vertex association
      double minVtxDist3D = 9999.;
      int indexVtx = -1;
      double vertexDistZ = 9999.;
 
      if (primaryVertices.isValid()) {
	edm::LogInfo("ElectronBlock") << "Total # Primary Vertices: " << primaryVertices->size();
        for (reco::VertexCollection::const_iterator v_it = primaryVertices->begin() ; 
                                                   v_it != primaryVertices->end() ; ++v_it) {
          double dist3D = sqrt(pow(it->gsfTrack()->dxy(v_it->position()),2) + pow(it->gsfTrack()->dz(v_it->position()),2));
          if (dist3D < minVtxDist3D) {
            minVtxDist3D = dist3D;
            indexVtx = int(std::distance(primaryVertices->begin(),v_it));
            vertexDistZ = it->gsfTrack()->dz(v_it->position());
          }
        }
      } 
      else {
	edm::LogError("ElectronBlock") << "Error! Can't get the product " << _vtxInputTag;
      }      
      double pfreliso = it->userFloat("pfLooseIsoPt04")/it->pt();

      electronB = new ((*cloneElectron)[fnElectron++]) vhtm::Electron();
      electronB->eta         = it->eta();
      electronB->phi         = it->phi();
      electronB->pt          = it->pt();
      bool hasGsfTrack = it->gsfTrack().isNonnull() ? true : false;
      electronB->hasGsfTrack = hasGsfTrack;
      electronB->trackPt     = (hasGsfTrack) ? it->gsfTrack()->pt() : -1;
      electronB->energy      = it->energy();
      electronB->caloEnergy  = it->caloEnergy();
      electronB->charge      = it->charge();
      electronB->nValidHits  = (hasGsfTrack) ? it->gsfTrack()->numberOfValidHits() : -1; 

      electronB->simpleEleId60cIso 
                             = it->electronID("simpleEleId60cIso");
      electronB->simpleEleId70cIso 
                             = it->electronID("simpleEleId70cIso");
      electronB->simpleEleId80cIso 
                             = it->electronID("simpleEleId80cIso");
      electronB->simpleEleId85cIso 
                             = it->electronID("simpleEleId85cIso");
      electronB->simpleEleId90cIso 
                             = it->electronID("simpleEleId90cIso");
      electronB->simpleEleId95cIso 
                             = it->electronID("simpleEleId95cIso");

      // ID variables
      electronB->hoe           = it->hadronicOverEm();
      electronB->eop           = it->eSuperClusterOverP(); 
      electronB->sigmaEtaEta   = it->sigmaEtaEta();
      electronB->sigmaIEtaIEta = it->sigmaIetaIeta();
      electronB->deltaPhiTrkSC = it->deltaPhiSuperClusterTrackAtVtx();
      electronB->deltaEtaTrkSC = it->deltaEtaSuperClusterTrackAtVtx();
      electronB->classif       = it->classification();
      electronB->e1x5overe5x5  = (it->e5x5() > 0) ? (it->e1x5()/it->e5x5()) : 0;
      electronB->e2x5overe5x5  = (it->e5x5() > 0) ? (it->e2x5Max()/it->e5x5()) : 0;

      // Iso variables
      electronB->isoEcal03 = it->dr03EcalRecHitSumEt();
      electronB->isoHcal03 = it->dr03HcalTowerSumEt();
      electronB->isoTrk03  = it->dr03TkSumPt();
      electronB->isoEcal04 = it->dr04EcalRecHitSumEt(); // ecalIso
      electronB->isoHcal04 = it->dr04HcalTowerSumEt(); // hcalIso
      electronB->isoTrk04  = it->dr04TkSumPt(); // trackIso
      electronB->isoRel03  = (it->dr03EcalRecHitSumEt()+it->dr03HcalTowerSumEt()+it->dr03TkSumPt())/it->pt();
      electronB->isoRel04  = (it->dr04EcalRecHitSumEt()+it->dr04HcalTowerSumEt()+it->dr04TkSumPt())/it->pt();

      // Conversion variables
      electronB->missingHits = (hasGsfTrack) ? it->gsfTrack()->trackerExpectedHitsInner().numberOfHits() : -1;
      electronB->dist_vec    = dist;
      electronB->dCotTheta   = dcot;

      // SC associated with electron
      electronB->scEn        = it->superCluster()->energy();
      electronB->scEta       = it->superCluster()->eta();
      electronB->scPhi       = it->superCluster()->phi();
      electronB->scET        = it->superCluster()->energy()/cosh(it->superCluster()->eta());
      electronB->scRawEnergy = it->superCluster()->rawEnergy();

      // Vertex association variables
      electronB->vtxDist3D = minVtxDist3D;
      electronB->vtxIndex  = indexVtx;
      electronB->vtxDistZ  = vertexDistZ;
      electronB->relIso    = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();
      electronB->pfRelIso  = pfreliso;

      // IP information
      electronB->dB  = it->dB(pat::Electron::PV2D);
      electronB->edB = it->edB(pat::Electron::PV2D);

      // Bremstrahlung information
      electronB->nBrems = it->numberOfBrems();
      electronB->fbrem  = it->fbrem();

      // MIT Electron ID
      EcalClusterLazyTools lazyTools(iEvent, iSetup, _ecalEBInputTag, _ecalEEInputTag);
      const pat::Electron& ele = *it;
      electronB->mva = fMVA->MVAValue(&ele,  lazyTools);
    }
  } 
  else {
    edm::LogError("ElectronBlock") << "Error! Can't get the product " << _electronInputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronBlock);
