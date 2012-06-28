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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Utilities/General/interface/FileInPath.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "VHTauTau/TreeMaker/plugins/ElectronBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

// Constructor
ElectronBlock::ElectronBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _bsInputTag(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
  _trkInputTag(iConfig.getParameter<edm::InputTag>("trackSrc")),
  _dcsInputTag(iConfig.getParameter<edm::InputTag>("dcsSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  _convInputTag(iConfig.getParameter<edm::InputTag>("convSrc")),
  _electronInputTag(iConfig.getParameter<edm::InputTag>("electronSrc")),
  _pfElectronInputTag(iConfig.getParameter<edm::InputTag>("pfElectronSrc")),
  _ecalEBInputTag(iConfig.getParameter<edm::InputTag>("ecalEBSrc")),
  _ecalEEInputTag(iConfig.getParameter<edm::InputTag>("ecalEESrc")),
  _rhoInputTag(iConfig.getParameter<edm::InputTag>("rhoSrc")),
  _pfInputTag(iConfig.getParameter<edm::InputTag>("pfSrc"))
{
  edm::FileInPath idCat1Weights = iConfig.getParameter<edm::FileInPath>("IdCat1Weights");
  edm::FileInPath idCat2Weights = iConfig.getParameter<edm::FileInPath>("IdCat2Weights");
  edm::FileInPath idCat3Weights = iConfig.getParameter<edm::FileInPath>("IdCat3Weights");
  edm::FileInPath idCat4Weights = iConfig.getParameter<edm::FileInPath>("IdCat4Weights");
  edm::FileInPath idCat5Weights = iConfig.getParameter<edm::FileInPath>("IdCat5Weights");
  edm::FileInPath idCat6Weights = iConfig.getParameter<edm::FileInPath>("IdCat6Weights");

  // Electron ID MVA
  std::vector<std::string> eleid_wtfiles;
  eleid_wtfiles.push_back(idCat1Weights.fullPath());
  eleid_wtfiles.push_back(idCat2Weights.fullPath());
  eleid_wtfiles.push_back(idCat3Weights.fullPath());
  eleid_wtfiles.push_back(idCat4Weights.fullPath());
  eleid_wtfiles.push_back(idCat5Weights.fullPath());
  eleid_wtfiles.push_back(idCat6Weights.fullPath());

  fElectronIdMVA = new EGammaMvaEleEstimator();
  fElectronIdMVA->initialize("BDT",
                        EGammaMvaEleEstimator::kNonTrig,
                        true, 
                        eleid_wtfiles);

  edm::FileInPath isoBarrelPt1Weights = iConfig.getParameter<edm::FileInPath>("IsoBarrelPt1Weights");
  edm::FileInPath isoEndcapPt1Weights = iConfig.getParameter<edm::FileInPath>("IsoEndcapPt1Weights");
  edm::FileInPath isoBarrelPt2Weights = iConfig.getParameter<edm::FileInPath>("IsoBarrelPt2Weights");
  edm::FileInPath isoEndcapPt2Weights = iConfig.getParameter<edm::FileInPath>("IsoEndcapPt2Weights");

  std::vector<std::string> eleiso_wtfiles;
  eleiso_wtfiles.push_back(isoBarrelPt1Weights.fullPath());
  eleiso_wtfiles.push_back(isoEndcapPt1Weights.fullPath());
  eleiso_wtfiles.push_back(isoBarrelPt2Weights.fullPath());
  eleiso_wtfiles.push_back(isoEndcapPt2Weights.fullPath());

  std::string target = iConfig.getParameter<std::string>("target");
  if (target == "2011Data") {
    target_ = ElectronEffectiveArea::kEleEAData2011;
  } else if (target == "2012Data") {
    target_ = ElectronEffectiveArea::kEleEAData2012;
  } else if (target == "Fall11MC") {
    target_ = ElectronEffectiveArea::kEleEAFall11MC;
  } else if (target == "Summer11MC") {
    target_ = ElectronEffectiveArea::kEleEASummer11MC;
  } 
  else {
    throw cms::Exception("UnknownTarget")
      << "Bad eff. area option for electrons: " << target
      << " options are: 2011Data, 2012Data, Fall11MC, Summer11MC" << std::endl;
  }
  fElectronIsoMVA = new EGammaMvaEleEstimator();
  fElectronIsoMVA->initialize("EleIso_BDTG_IsoRings",
			      EGammaMvaEleEstimator::kIsoRings,
			      true,
			      eleiso_wtfiles);
}
ElectronBlock::~ElectronBlock() { 
  delete fElectronIdMVA; 
  delete fElectronIsoMVA; 
}
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

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel(_bsInputTag, beamSpot);

  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByLabel(_dcsInputTag, dcsHandle);

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel(_convInputTag, hConversions);

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(_electronInputTag, electrons);

  edm::Handle<reco::PFCandidateCollection> pfElectrons;
  iEvent.getByLabel(_pfElectronInputTag, pfElectrons);

  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  TransientTrackBuilder ttrackBuilder = *(builder.product());

  edm::Handle<double> hRho;
  iEvent.getByLabel(_rhoInputTag, hRho);
  double Rho = *hRho;

  edm::Handle<reco::PFCandidateCollection> hPfCandProduct;
  iEvent.getByLabel(_pfInputTag, hPfCandProduct);
  const reco::PFCandidateCollection &inPfCands = *(hPfCandProduct.product());

  // Just leave these blank.
  reco::GsfElectronCollection identifiedElectrons;
  reco::MuonCollection identifiedMuons;

  double evt_bField = 3.8;
  // need the magnetic field
  //
  // if isRealData then derive bfield using the
  // magnet current from DcsStatus
  // otherwise take it from the IdealMagneticFieldRecord
  if (iEvent.isRealData()) {
    if (dcsHandle.isValid()) {
      edm::LogInfo("ElectronBlock") << "Success >> obtained DcsStatusCollection for label:" 
                                    << _dcsInputTag;
      // scale factor = 3.801/18166.0 which are
      // average values taken over a stable two-week period
      double currentToBFieldScaleFactor = 2.09237036221512717e-04;
      double current = (*dcsHandle)[0].magnetCurrent();
      evt_bField = current*currentToBFieldScaleFactor;
    } 
    else {
      edm::LogError("ElectronBlock") << "Error >> Failed to get DcsStatusCollection for label: " 
                                     << _dcsInputTag;
    }
  } 
  else {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    if (magneticField.isValid()) {
      edm::LogInfo("ElectronBlock") << "Success >> obtained IdealMagneticFieldRecord!";
      evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
    } 
    else {
      edm::LogError("ElectronBlock") << "Error >> Failed to get IdealMagneticFieldRecord";
    }
  }

  if (electrons.isValid()) {
    edm::LogInfo("ElectronBlock") << "Total # PAT Electrons: " << electrons->size();
    for (std::vector<pat::Electron>::const_iterator it  = electrons->begin(); 
                                                    it != electrons->end(); ++it) {
      if (fnElectron == kMaxElectron) {
	edm::LogInfo("ElectronBlock") << "Too many PAT Electrons, fnElectron = " 
                                      << fnElectron; 
	break;
      }
      // if electron is not ECAL driven, continue
      if (!it->ecalDrivenSeed()) continue;

      bool hasGsfTrack  = it->gsfTrack().isNonnull() ? true : false;
      reco::GsfTrackRef tk = it->gsfTrack();

      electronB = new ((*cloneElectron)[fnElectron++]) vhtm::Electron();
      electronB->eta         = it->eta();
      electronB->phi         = it->phi();
      electronB->pt          = it->pt();
      electronB->hasGsfTrack = hasGsfTrack;
      electronB->energy      = it->energy();
      electronB->caloEnergy  = it->ecalEnergy();
      electronB->caloEnergyError = it->ecalEnergyError();
      electronB->charge      = it->charge();
  
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

      if (hasGsfTrack) {
        electronB->trackPt      = tk->pt();
        electronB->trackPtError = tk->ptError();

	// Hit pattern
	const reco::HitPattern& hitp = tk->hitPattern();
	electronB->pixHits = hitp.numberOfValidPixelHits();
	electronB->trkHits = hitp.numberOfValidTrackerHits();

        electronB->nValidHits   = tk->numberOfValidHits(); 
        electronB->missingHits  = tk->trackerExpectedHitsInner().numberOfHits();

        electronB->trkD0        = tk->d0();
        electronB->trkD0Error   = tk->d0Error();
      }
      // ID variables
      electronB->hoe           = it->hcalOverEcal();
      electronB->hoeDepth1     = it->hcalDepth1OverEcal();
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
      electronB->isoRel03  = (it->dr03EcalRecHitSumEt()
                            + it->dr03HcalTowerSumEt()
                            + it->dr03TkSumPt())/it->pt();
      electronB->isoRel04  = (it->dr04EcalRecHitSumEt()
                            + it->dr04HcalTowerSumEt()
                            + it->dr04TkSumPt())/it->pt();

      // Conversion variables
      ConversionFinder convFinder;
      double dist = -9999.;
      double dcot = -9999.;
      if (tracks.isValid()) {
	edm::LogInfo("ElectronBlock") << "Success >> obtained TrackCollection for label: " 
                                      << _trkInputTag;
        ConversionInfo convInfo = convFinder.getConversionInfo(*it, tracks, evt_bField);
        dist = convInfo.dist();
        dcot = convInfo.dcot();
      } 
      else {
	edm::LogError("ElectronBlock") << "Error >> Failed to get TrackCollection for label: " 
                                       << _trkInputTag;
      }
      electronB->dist_vec  = dist;
      electronB->dCotTheta = dcot;

      const reco::GsfElectron* aGsf = static_cast<const reco::GsfElectron*>(&(*it));
      bool hasMatchedConv = ConversionTools::hasMatchedConversion(*aGsf, 
                                                                  hConversions, 
                                                                  beamSpot->position(),
                                                                  true, 2.0, 1e-06,0);
      electronB->hasMatchedConv = hasMatchedConv;

      // SC associated with electron
      electronB->scEn        = it->superCluster()->energy();
      electronB->scEta       = it->superCluster()->eta();
      electronB->scPhi       = it->superCluster()->phi();
      electronB->scET        = it->superCluster()->energy()/cosh(it->superCluster()->eta());
      electronB->scRawEnergy = it->superCluster()->rawEnergy();

      // Vertex association
      double minVtxDist3D = 9999.;
      int indexVtx = -1;
      double vertexDistZ = 9999.;
      if (hasGsfTrack) {
        if (primaryVertices.isValid()) {
	  edm::LogInfo("ElectronBlock") << "Total # Primary Vertices: " << primaryVertices->size();
          for (reco::VertexCollection::const_iterator vit  = primaryVertices->begin(); 
                                                      vit != primaryVertices->end(); ++vit) {
            double dxy = tk->dxy(vit->position());
            double dz  = tk->dz(vit->position());
            double dist3D = std::sqrt(pow(dxy, 2) + pow(dz, 2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = int(std::distance(primaryVertices->begin(), vit));
              vertexDistZ = dz;
            }
          }
        } 
        else {
	  edm::LogError("ElectronBlock") << "Error >> Failed to get VertexCollection for label: " 
                                         << _vtxInputTag;
        }      
      }
      // Vertex association variables
      electronB->vtxDist3D = minVtxDist3D;
      electronB->vtxIndex  = indexVtx;
      electronB->vtxDistZ  = vertexDistZ;

      // UW prescription for pf based isolation
      double v1 = it->photonIso() + it->neutralHadronIso() - 0.5*it->userIso(0);
      double vmax = (v1 > 0) ? v1 : 0;
      double UWpfreliso = (it->chargedHadronIso() + vmax)/it->pt();

      electronB->relIso   = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();
      electronB->pfRelIso = UWpfreliso;

      // PFlow isolation information
      electronB->chargedHadronIso = it->pfIsolationVariables().chargedHadronIso;
      electronB->neutralHadronIso = it->pfIsolationVariables().neutralHadronIso;
      electronB->photonIso        = it->pfIsolationVariables().photonIso;
  
      // IP information
      electronB->dB  = it->dB(pat::Electron::PV2D);
      electronB->edB = it->edB(pat::Electron::PV2D);

      electronB->dB3d  = it->dB(pat::Electron::PV3D);
      electronB->edB3d = it->edB(pat::Electron::PV3D);

      // Bremstrahlung information
      electronB->nBrems = it->numberOfBrems();
      electronB->fbrem  = it->fbrem();

      // MIT MVA Electron ID
      EcalClusterLazyTools lazyTools(iEvent, iSetup, _ecalEBInputTag, _ecalEEInputTag);
      //const pat::Electron& ele = *it;
      double idMVA = (primaryVertices->size()) 
	? fElectronIdMVA->mvaValue(*aGsf, primaryVertices->at(0), ttrackBuilder, lazyTools, false) // &ele
        : -1;
      electronB->idMVA = idMVA;

      // Isolation MVA
      double isomva = (primaryVertices->size()) 
        ? fElectronIsoMVA->mvaValue(*aGsf, primaryVertices->at(0),
				    inPfCands, Rho,
				    target_,
				    identifiedElectrons, identifiedMuons)
        : -1;
      electronB->isoMVA = isomva;

      // Fiducial flag 
      int fidFlag = 0;
      if (it->isEB())        fidFlag |= (1 << 0);
      if (it->isEE())        fidFlag |= (1 << 1);
      if (it->isEBEtaGap())  fidFlag |= (1 << 2);
      if (it->isEBPhiGap())  fidFlag |= (1 << 3);
      if (it->isEERingGap()) fidFlag |= (1 << 4);
      if (it->isEEDeeGap())  fidFlag |= (1 << 5);
      if (it->isEBEEGap())   fidFlag |= (1 << 6);
      electronB->fidFlag = fidFlag;

      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      electronB->vx = vertex.x();             
      electronB->vy = vertex.y();             
      electronB->vz = vertex.z();             

      // Iso deposit and PF Isolaiton
      fillIsoDeposit(*it, electronB);
    }
  } 
  else {
    edm::LogError("ElectronBlock") << "Error >> Failed to get pat::Electron Collection for label: " 
                                   << _electronInputTag;
  }
}
void ElectronBlock::fillIsoDeposit(const pat::Electron& ele, vhtm::Electron* electronB) {
  double pt = ele.pt();
  double eta = ele.eta();
  double phi = ele.phi();

  reco::isodeposit::AbsVetos v10Charged;
  reco::isodeposit::AbsVetos v10Neutral;  
  reco::isodeposit::AbsVetos v10Photons;
  reco::isodeposit::AbsVetos v11Charged; 
  reco::isodeposit::AbsVetos v11Neutral;  
  reco::isodeposit::AbsVetos v11Photons;
  reco::isodeposit::AbsVetos v11EECharged; 
  reco::isodeposit::AbsVetos v11EENeutral;  
  reco::isodeposit::AbsVetos v11EEPhotons;

  v10Charged.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi),0.01));
  v10Charged.push_back(new reco::isodeposit::ThresholdVeto(0.5));

  v10Neutral.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta,phi),0.08));
  v10Neutral.push_back(new reco::isodeposit::ThresholdVeto(1.0));

  v10Photons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi),0.05));
  v10Photons.push_back(new reco::isodeposit::ThresholdVeto(1.0));

  v11Charged.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi),0.03));
  v11Charged.push_back(new reco::isodeposit::ThresholdVeto(0.5));

  v11Neutral.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi),0.08));
  v11Neutral.push_back(new reco::isodeposit::ThresholdVeto(0.5));

  v11Photons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi),0.05));
  v11Photons.push_back(new reco::isodeposit::ThresholdVeto(0.5));

  v11EECharged.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi),0.03));
  v11EECharged.push_back(new reco::isodeposit::ThresholdVeto(0.5));

  v11EENeutral.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi),0.08));
  v11EENeutral.push_back(new reco::isodeposit::ThresholdVeto(0.5));

  v11EEPhotons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi),0.05));
  v11EEPhotons.push_back(new reco::isodeposit::ThresholdVeto(0.5));

  float chIso03v1   = ele.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, v10Charged).first;
  float nhIso03v1   = ele.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, v10Neutral).first;
  float phIso03v1   = ele.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, v10Photons).first;
  float nhIsoPU03v1 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v10Neutral).first;
  float phIsoPU03v1 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v10Photons).first;

  float chIso04v1   = ele.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, v10Charged).first;
  float nhIso04v1   = ele.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, v10Neutral).first;
  float phIso04v1   = ele.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, v10Photons).first;
  float nhIsoPU04v1 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v10Neutral).first;
  float phIsoPU04v1 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v10Photons).first;

  float chIso03v2   = ele.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, v11Charged).first;
  float nhIso03v2   = ele.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, v11Neutral).first;
  float phIso03v2   = ele.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, v11Photons).first;
  float nhIsoPU03v2 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v11Neutral).first;
  float phIsoPU03v2 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v11Photons).first;

  float chIso04v2   = ele.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, v11Charged).first;
  float nhIso04v2   = ele.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, v11Neutral).first;
  float phIso04v2   = ele.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, v11Photons).first;
  float nhIsoPU04v2 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v11Neutral).first;
  float phIsoPU04v2 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v11Photons).first;

  float chIso03EEv2   = ele.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, v11EECharged).first;
  float nhIso03EEv2   = ele.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, v11EENeutral).first;
  float phIso03EEv2   = ele.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, v11EEPhotons).first;
  float nhIsoPU03EEv2 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v11EENeutral).first;
  float phIsoPU03EEv2 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v11EEPhotons).first;

  float chIso04EEv2   = ele.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, v11EECharged).first;
  float nhIso04EEv2   = ele.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, v11EENeutral).first;
  float phIso04EEv2   = ele.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, v11EEPhotons).first;
  float nhIsoPU04EEv2 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v11EENeutral).first;
  float phIsoPU04EEv2 = ele.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v11EEPhotons).first;

  int flagEB = (ele.isEB() ? 1 : 0);
  int flagEE = (ele.isEE() ? 1 : 0); 

  electronB->pfRelIso03v1   = (chIso03v1 + nhIso03v1 + phIso03v1)/pt;
  electronB->pfRelIsoDB03v1 = (chIso03v1 + std::max(nhIso03v1 + phIso03v1 - 0.5 * 0.5 * (nhIsoPU03v1 + phIsoPU03v1), 0.0))/pt;
  electronB->pfRelIso03v2   = (chIso03v2 + nhIso03v2 + phIso03v2)/pt;
  electronB->pfRelIsoDB03v2 = (chIso03v2 + std::max(nhIso03v2 + phIso03v2 - 0.5 * 0.5 * (nhIsoPU03v2 + phIsoPU03v2), 0.0))/pt;
  electronB->pfRelIsoDB03v3 =
        flagEB * (chIso03v2   + std::max(nhIso03v2   + phIso03v2   - 0.5 * 0.5 * (nhIsoPU03v2 + phIsoPU03v2), 0.0))/pt
      + flagEE * (chIso03EEv2 + std::max(nhIso03EEv2 + phIso03EEv2 - 0.5 * 0.5 * (nhIsoPU03EEv2+phIsoPU03EEv2),0.0))/pt;
    
  electronB->pfRelIso04v1   = (chIso04v1 + nhIso04v1 + phIso04v1)/pt;
  electronB->pfRelIsoDB04v1 = (chIso04v1 + std::max(nhIso04v1 + phIso04v1 - 0.5 * 0.5 * (nhIsoPU04v1 + phIsoPU04v1), 0.0))/pt;
  electronB->pfRelIso04v2   = (chIso04v2 + nhIso04v2 + phIso04v2)/pt;
  electronB->pfRelIsoDB04v2 = (chIso04v2 + std::max(nhIso04v2 + phIso04v2 - 0.5 * 0.5 * (nhIsoPU04v2 + phIsoPU04v2), 0.0))/pt;
  electronB->pfRelIsoDB04v3 =
        flagEB * (chIso04v2   + std::max(nhIso04v2   + phIso04v2   - 0.5 * 0.5 * (nhIsoPU04v2   + phIsoPU04v2), 0.0))/pt
      + flagEE * (chIso04EEv2 + std::max(nhIso04EEv2 + phIso04EEv2 - 0.5 * 0.5 * (nhIsoPU04EEv2 + phIsoPU04EEv2), 0.0))/pt;

  // cleaning
  for (std::vector<reco::isodeposit::AbsVeto*>::const_iterator it  = v10Charged.begin();
                                                               it != v10Charged.end(); ++it) {
    delete (*it);
  }
  for (std::vector<reco::isodeposit::AbsVeto*>::const_iterator it  = v10Neutral.begin();
                                                               it != v10Neutral.end(); ++it) {
    delete (*it);
  }
  for (std::vector<reco::isodeposit::AbsVeto*>::const_iterator it  = v10Photons.begin();
                                                               it != v10Photons.end(); ++it) {
    delete (*it);
  }
  for (std::vector<reco::isodeposit::AbsVeto*>::const_iterator it  = v11Charged.begin();
                                                               it != v11Charged.end(); ++it) {
    delete (*it);
  }
  for (std::vector<reco::isodeposit::AbsVeto*>::const_iterator it  = v11Neutral.begin();
                                                               it != v11Neutral.end(); ++it) {
    delete (*it);
  }
  for (std::vector<reco::isodeposit::AbsVeto*>::const_iterator it  = v11Photons.begin();
                                                               it != v11Photons.end(); ++it) {
    delete (*it);
  }
  for (std::vector<reco::isodeposit::AbsVeto*>::const_iterator it  = v11EECharged.begin();
                                                               it != v11EECharged.end(); ++it) {
    delete (*it);
  }
  for (std::vector<reco::isodeposit::AbsVeto*>::const_iterator it  = v11EENeutral.begin();
                                                               it != v11EENeutral.end(); ++it) {
    delete (*it);
  }
  for (std::vector<reco::isodeposit::AbsVeto*>::const_iterator it  = v11EEPhotons.begin();
                                                               it != v11EEPhotons.end(); ++it) {
    delete (*it);
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronBlock);
