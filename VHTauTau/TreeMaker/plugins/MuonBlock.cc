#include <iostream>

#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"

#include "VHTauTau/TreeMaker/plugins/MuonBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

MuonBlock::MuonBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _muonInputTag(iConfig.getParameter<edm::InputTag>("muonSrc")),
  _pfMuonInputTag(iConfig.getParameter<edm::InputTag>("pfMuonSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  _beamSpotInputTag(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
  _beamSpotCorr(iConfig.getParameter<bool>("beamSpotCorr")),
  _muonID(iConfig.getParameter<std::string>("muonID"))
{
  // ID
  edm::FileInPath idBarrel1Weights = iConfig.getParameter<edm::FileInPath>("IdBarrel1Weights");
  edm::FileInPath idEndcap1Weights = iConfig.getParameter<edm::FileInPath>("IdEndcap1Weights");
  edm::FileInPath idBarrel2Weights = iConfig.getParameter<edm::FileInPath>("IdBarrel2Weights");
  edm::FileInPath idEndcap2Weights = iConfig.getParameter<edm::FileInPath>("IdEndcap2Weights");
  edm::FileInPath idTrackerWeights = iConfig.getParameter<edm::FileInPath>("IdTrackerWeights");
  edm::FileInPath idGlobalWeights  = iConfig.getParameter<edm::FileInPath>("IdGlobalWeights");

  std::vector<std::string> muonid_wtfiles;
  muonid_wtfiles.push_back(idBarrel1Weights.fullPath());
  muonid_wtfiles.push_back(idBarrel2Weights.fullPath());
  muonid_wtfiles.push_back(idEndcap1Weights.fullPath());
  muonid_wtfiles.push_back(idEndcap2Weights.fullPath());
  muonid_wtfiles.push_back(idTrackerWeights.fullPath());
  muonid_wtfiles.push_back(idGlobalWeights.fullPath());

  fMuonIdMVA = new MuonMVAEstimator();
  fMuonIdMVA->initialize("MuonID_BDTG",
                         MuonMVAEstimator::kID,
                         true,
                         muonid_wtfiles);

  // Rings isolation
  edm::FileInPath isoBarrel1Weights = iConfig.getParameter<edm::FileInPath>("IsoBarrel1Weights");
  edm::FileInPath isoEndcap1Weights = iConfig.getParameter<edm::FileInPath>("IsoEndcap1Weights");
  edm::FileInPath isoBarrel2Weights = iConfig.getParameter<edm::FileInPath>("IsoBarrel2Weights");
  edm::FileInPath isoEndcap2Weights = iConfig.getParameter<edm::FileInPath>("IsoEndcap2Weights");
  edm::FileInPath isoTrackerWeights = iConfig.getParameter<edm::FileInPath>("IsoTrackerWeights");
  edm::FileInPath isoGlobalWeights  = iConfig.getParameter<edm::FileInPath>("IsoGlobalWeights");

  std::vector<std::string> muoniso_wtfiles;
  muoniso_wtfiles.push_back(isoBarrel1Weights.fullPath());
  muoniso_wtfiles.push_back(isoEndcap1Weights.fullPath());
  muoniso_wtfiles.push_back(isoBarrel2Weights.fullPath());
  muoniso_wtfiles.push_back(isoEndcap2Weights.fullPath());
  muoniso_wtfiles.push_back(isoTrackerWeights.fullPath());
  muoniso_wtfiles.push_back(isoGlobalWeights.fullPath());

  fMuonIsoMVA = new MuonMVAEstimator();
  fMuonIsoMVA->initialize("MuonIso_BDTG_IsoRings",
                          MuonMVAEstimator::kIsoRings,
                          true,
                          muoniso_wtfiles);

  // Radial + Rings ISO
  edm::FileInPath isoRingsRadBarrel1Weights = iConfig.getParameter<edm::FileInPath>("IsoRingsRadBarrel1Weights");
  edm::FileInPath isoRingsRadEndcap1Weights = iConfig.getParameter<edm::FileInPath>("IsoRingsRadEndcap1Weights");
  edm::FileInPath isoRingsRadBarrel2Weights = iConfig.getParameter<edm::FileInPath>("IsoRingsRadBarrel2Weights");
  edm::FileInPath isoRingsRadEndcap2Weights = iConfig.getParameter<edm::FileInPath>("IsoRingsRadEndcap2Weights");
  edm::FileInPath isoRingsRadTrackerWeights = iConfig.getParameter<edm::FileInPath>("IsoRingsRadTrackerWeights");
  edm::FileInPath isoRingsRadGlobalWeights  = iConfig.getParameter<edm::FileInPath>("IsoRingsRadGlobalWeights");

  std::vector<std::string> muonisoRingsRad_wtfiles;
  muonisoRingsRad_wtfiles.push_back(isoRingsRadBarrel1Weights.fullPath());
  muonisoRingsRad_wtfiles.push_back(isoRingsRadEndcap1Weights.fullPath());
  muonisoRingsRad_wtfiles.push_back(isoRingsRadBarrel2Weights.fullPath());
  muonisoRingsRad_wtfiles.push_back(isoRingsRadEndcap2Weights.fullPath());
  muonisoRingsRad_wtfiles.push_back(isoRingsRadTrackerWeights.fullPath());
  muonisoRingsRad_wtfiles.push_back(isoRingsRadGlobalWeights.fullPath());

  fMuonIsoRingsRadMVA = new MuonMVAEstimator();
  fMuonIsoRingsRadMVA->initialize("MuonIso_BDTG_IsoRingsRad",
                                  MuonMVAEstimator::kIsoRingsRadial,
                                  true,
                                  muonisoRingsRad_wtfiles);

  std::string target = iConfig.getParameter<std::string>("target");
  if (target == "2011Data") {
    target_ = MuonEffectiveArea::kMuEAData2011;
  } else if (target == "2012Data") {
    target_ = MuonEffectiveArea::kMuEAData2012;
  } else if (target == "Fall11MC") {
    target_ = MuonEffectiveArea::kMuEAFall11MC;
  } else if (target == "Summer11MC") {
    target_ = MuonEffectiveArea::kMuEASummer11MC;
  } else {
    throw cms::Exception("UnknownTarget")
      << "Bad eff. area option for muons: " << target
      << " options are: 2011Data, 2012Data, Fall11MC, Summer11MC" << std::endl;
  }
}
MuonBlock::~MuonBlock() {
  delete fMuonIdMVA;
  delete fMuonIsoMVA;
  delete fMuonIsoRingsRadMVA;
}
void MuonBlock::beginJob() 
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneMuon = new TClonesArray("vhtm::Muon");
  tree->Branch("Muon", &cloneMuon, 32000, 2);
  tree->Branch("nMuon", &fnMuon,  "fnMuon/I");
}
void MuonBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneMuon->Clear();
  fnMuon = 0;

  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(_muonInputTag, muons);

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel(_beamSpotInputTag, beamSpot);

  edm::Handle<reco::PFCandidateCollection> pfMuons;
  iEvent.getByLabel(_pfMuonInputTag, pfMuons);

  edm::Handle<double> hRho;
  iEvent.getByLabel(_rhoInputTag, hRho);
  double Rho = *hRho;

  edm::Handle<reco::PFCandidateCollection> hPfCandProduct;
  iEvent.getByLabel(_pfInputTag, hPfCandProduct);
  const reco::PFCandidateCollection &inPfCands = *(hPfCandProduct.product());

  // Just leave these blank.
  reco::GsfElectronCollection identifiedElectrons;
  reco::MuonCollection identifiedMuons;

  if (muons.isValid()) {
    edm::LogInfo("MuonBlock") << "Total # Muons: " << muons->size();
    for (std::vector<pat::Muon>::const_iterator it  = muons->begin(); 
                                                it != muons->end(); ++it) {
      // consider only global muons
      if (!it->isGlobalMuon()) continue;

      if (fnMuon == kMaxMuon) {
	edm::LogInfo("MuonBlock") << "Too many PAT Muons, fnMuon = " << fnMuon; 
	break;
      }

      reco::TrackRef tk  = it->innerTrack();  // tracker segment only
      reco::TrackRef gtk = it->globalTrack(); 

      muonB = new ((*cloneMuon)[fnMuon++]) vhtm::Muon();
      muonB->isTrackerMuon = (it->isTrackerMuon()) ? true : false;

      muonB->eta        = it->eta();
      muonB->phi        = it->phi();
      muonB->pt         = it->pt();
      muonB->ptError    = tk->ptError();
      muonB->p          = it->p();
      muonB->energy     = it->energy();
      muonB->charge     = it->charge();

      double trkd0 = tk->d0();
      double trkdz = tk->dz();
      if (_beamSpotCorr) {
        if (beamSpot.isValid()) {
          trkd0 = -(tk->dxy(beamSpot->position()));
          trkdz = tk->dz(beamSpot->position());
        }
        else
          edm::LogError("MuonsBlock") << "Error >> Failed to get BeamSpot for label: "
                                      << _beamSpotInputTag;
      }
      muonB->trkD0      = trkd0;
      muonB->trkD0Error = tk->d0Error();
      muonB->trkDz      = trkdz;
      muonB->trkDzError = tk->dzError();
      muonB->globalChi2 = it->normChi2();
      muonB->passID     = (it->muonID(_muonID)) ? true : false;

      // Vertex association
      double minVtxDist3D = 9999.;
         int indexVtx = -1;
      double vertexDistZ = 9999.;
      if (primaryVertices.isValid()) {
	edm::LogInfo("MuonBlock") << "Total # Primary Vertices: " << primaryVertices->size();

        for (reco::VertexCollection::const_iterator vit  = primaryVertices->begin(); 
                                                    vit != primaryVertices->end(); ++vit) {
          double dxy = tk->dxy(vit->position());
          double dz  = tk->dz(vit->position());
          double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
          if (dist3D < minVtxDist3D) {
            minVtxDist3D = dist3D;
            indexVtx = int(std::distance(primaryVertices->begin(), vit));
            vertexDistZ = dz;
          }
        }
      } 
      else {
	edm::LogError("MuonBlock") << "Error >> Failed to get VertexCollection for label: " 
                                   << _vtxInputTag;
      }
      muonB->vtxDist3D = minVtxDist3D;
      muonB->vtxIndex  = indexVtx;
      muonB->vtxDistZ  = vertexDistZ;

      // Hit pattern
      const reco::HitPattern& hitp = gtk->hitPattern();  // innerTrack will not provide Muon Hits 
      muonB->pixHits = hitp.numberOfValidPixelHits();
      muonB->trkHits = hitp.numberOfValidTrackerHits();
      muonB->muoHits = hitp.numberOfValidMuonHits();
      muonB->matches = it->numberOfMatches();

      // Isolation
      muonB->trkIso  = it->trackIso();
      muonB->ecalIso = it->ecalIso();
      muonB->hcalIso = it->hcalIso();
      muonB->hoIso   = it->isolationR03().hoEt;
      double reliso  = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();
      muonB->relIso  = reliso;

      // UW prescription for pf based isolation
      double v1 = it->photonIso() + it->neutralHadronIso() - 0.5*it->userIso(0);
      double vmax = (v1 > 0) ? v1 : 0;
      double UWpfreliso = (it->chargedHadronIso() + vmax)/it->pt();
      muonB->pfRelIso = UWpfreliso;

      // IP information
      muonB->dB  = it->dB(pat::Muon::PV2D);
      muonB->edB = it->edB(pat::Muon::PV2D);

      muonB->dB3d  = it->dB(pat::Muon::PV3D);
      muonB->edB3d = it->edB(pat::Muon::PV3D);

      // UW recommendation
      muonB->isGlobalMuonPromptTight = muon::isGoodMuon(*it, muon::GlobalMuonPromptTight);
      muonB->isAllArbitrated         = muon::isGoodMuon(*it, muon::AllArbitrated);
      muonB->nChambers               = it->numberOfChambers();
      muonB->nMatches                = it->numberOfMatches();
      muonB->nMatchedStations        = it->numberOfMatchedStations();
      muonB->stationMask             = it->stationMask();
      muonB->stationGapMaskDistance  = it->stationGapMaskDistance();
      muonB->stationGapMaskPull      = it->stationGapMaskPull();

      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      muonB->vx = vertex.x();             
      muonB->vy = vertex.y();             
      muonB->vz = vertex.z();             

      // Various MVA
      const pat::Muon& muon = *it;
      double idmva = fMuonIdMVA->mvaValue(muon, primaryVertices->at(0),
                                          inPfCands, Rho,
                                          target_,
                                          identifiedElectrons, identifiedMuons);
      
      double isomva = fMuonIsoMVA->mvaValue(muon, primaryVertices->at(0),
                                            inPfCands, Rho,
                                            target_,
                                            identifiedElectrons, identifiedMuons);

      double isoringsradmva 
         = fMuonIsoRingsRadMVA->mvaValue(muon, primaryVertices->at(0),
                                         inPfCands, Rho,
                                         target_,
                                         identifiedElectrons, identifiedMuons);
      muonB->idMVA  = idmva;
      muonB->isoMVA = isomva;
      muonB->isoRingsRadMVA = isoringsradmva;

      // Iso deposit and PF Isolaiton
      fillIsoDeposit(*it, muonB);
    }
  } 
  else {
    edm::LogError("MuonBlock") << "Error >> Failed to get pat::Muon collection for label: " 
                               << _muonInputTag;
  }
}
void MuonBlock::fillIsoDeposit(const pat::Muon& muo, vhtm::Muon* muonB) {
  double pt = muo.pt();
  double eta = muo.eta();
  double phi = muo.phi();

  reco::isodeposit::AbsVetos v10Charged;
  reco::isodeposit::AbsVetos v10Neutral;  
  reco::isodeposit::AbsVetos v10Photons;
  reco::isodeposit::AbsVetos v11Charged; 
  reco::isodeposit::AbsVetos v11Neutral;  
  reco::isodeposit::AbsVetos v11Photons;

  v10Charged.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.01));
  v10Charged.push_back(new reco::isodeposit::ThresholdVeto(0.5));
  v10Neutral.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.08));
  v10Neutral.push_back(new reco::isodeposit::ThresholdVeto(1.0));
  v10Photons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.05));
  v10Photons.push_back(new reco::isodeposit::ThresholdVeto(1.0));
    
  v11Charged.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.01));
  v11Charged.push_back(new reco::isodeposit::ThresholdVeto(0.5));
  v11Neutral.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.01));
  v11Neutral.push_back(new reco::isodeposit::ThresholdVeto(0.5));
  v11Photons.push_back(new reco::isodeposit::ConeVeto(reco::isodeposit::Direction(eta, phi), 0.01));
  v11Photons.push_back(new reco::isodeposit::ThresholdVeto(0.5));

  float chIso03v1   = muo.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, v10Charged).first;
  float nhIso03v1   = muo.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, v10Neutral).first;
  float phIso03v1   = muo.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, v10Photons).first;
  float nhIsoPU03v1 = muo.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v10Neutral).first;
  float phIsoPU03v1 = muo.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v10Photons).first;

  float chIso04v1   = muo.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, v10Charged).first;
  float nhIso04v1   = muo.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, v10Neutral).first;
  float phIso04v1   = muo.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, v10Photons).first;
  float nhIsoPU04v1 = muo.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v10Neutral).first;
  float phIsoPU04v1 = muo.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v10Photons).first;

  float chIso03v2   = muo.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, v11Charged).first;
  float nhIso03v2   = muo.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, v11Neutral).first;
  float phIso03v2   = muo.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, v11Photons).first;
  float nhIsoPU03v2 = muo.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v11Neutral).first;
  float phIsoPU03v2 = muo.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, v11Photons).first;

  float chIso04v2   = muo.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, v11Charged).first;
  float nhIso04v2   = muo.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, v11Neutral).first;
  float phIso04v2   = muo.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, v11Photons).first;
  float nhIsoPU04v2 = muo.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v11Neutral).first;
  float phIsoPU04v2 = muo.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, v11Photons).first;

  muonB->pfRelIso03v1   = (chIso03v1 + nhIso03v1 + phIso03v1)/pt;
  muonB->pfRelIso03v2   = (chIso03v2 + nhIso03v2 + phIso03v2)/pt;
  muonB->pfRelIsoDB03v1 = (chIso03v1 + std::max(nhIso03v1 + phIso03v1 - 0.5 * 0.5 * (nhIsoPU03v1 + phIsoPU03v1), 0.0))/pt;
  muonB->pfRelIsoDB03v2 = (chIso03v2 + std::max(nhIso03v2 + phIso03v2 - 0.5 * 0.5 * (nhIsoPU03v2 + phIsoPU03v2), 0.0))/pt;

  muonB->pfRelIso04v1   = (chIso04v1 + nhIso04v1 + phIso04v1)/pt;
  muonB->pfRelIso04v2   = (chIso04v2 + nhIso04v2 + phIso04v2)/pt;
  muonB->pfRelIsoDB04v1 = (chIso04v1 + std::max(nhIso04v1 + phIso04v1 - 0.5 * 0.5 * (nhIsoPU04v1 + phIsoPU04v1), 0.0))/pt;
  muonB->pfRelIsoDB04v2 = (chIso04v2 + std::max(nhIso04v2 + phIso04v2 - 0.5 * 0.5 * (nhIsoPU04v2 + phIsoPU04v2), 0.0))/pt;

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
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonBlock);
