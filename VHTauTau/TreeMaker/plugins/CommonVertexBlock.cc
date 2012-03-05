#include <iostream>
#include <vector>
#include <algorithm>

#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "VHTauTau/TreeMaker/plugins/CommonVertexBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

namespace {
  template <class T>
  class PtComparator {
  public:
    bool operator()(const T& a, const T& b) const {
      return a.obj()->pt() > b.obj()->pt();
    }
  };

  reco::TransientTrack getTracks(const reco::Candidate* cand, const TransientTrackBuilder* builder) {
    const pat::Muon* muon = dynamic_cast<const pat::Muon*>(cand);
    const pat::Electron* electron = dynamic_cast<const pat::Electron*>(cand);
    const pat::Tau* tau = dynamic_cast<const pat::Tau*>(cand);
    if (muon) {
      if (muon->innerTrack().isNonnull())
        return (builder->build(muon->innerTrack()));
    } 
    else if (electron) {
      if (electron->gsfTrack().isNonnull())
        return (builder->build(electron->gsfTrack()));
    } 
    else if (tau) {
      if (tau->signalPFChargedHadrCands()[0]->trackRef().isNonnull())
        return (builder->build(tau->signalPFChargedHadrCands()[0]->trackRef()));
      else
        if (tau->signalPFChargedHadrCands()[0]->gsfTrackRef().isNonnull())
          return (builder->build(tau->signalPFChargedHadrCands()[0]->gsfTrackRef()));
    }
    return reco::TransientTrack();
  }
}

CommonVertexBlock::CommonVertexBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _patMuonSrc(iConfig.getParameter<edm::InputTag>("patMuonSrc")),
  _patElectronSrc(iConfig.getParameter<edm::InputTag>("patElectronSrc")),
  _patTauSrc(iConfig.getParameter<edm::InputTag>("patTauSrc")),
  _minPtMuon(iConfig.getParameter<double>("minPtMuon")),
  _maxEtaMuon(iConfig.getParameter<double>("maxEtaMuon")),
  _maxChi2Muon(iConfig.getParameter<double>("maxChi2Muon")),
  _minTrkHitsMuon(iConfig.getParameter<double>("minTrkHitsMuon")),
  _minPtElectron(iConfig.getParameter<double>("minPtElectron")),
  _maxEtaElectron(iConfig.getParameter<double>("maxEtaElectron")),
  _minPtTau(iConfig.getParameter<double>("minPtTau"))
{}
void CommonVertexBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneCommonVertex = new TClonesArray("vhtm::CommonVertex");
  tree->Branch("CommonVertex", &cloneCommonVertex, 32000, 2);
  tree->Branch("nCommonVertex", &fnCommonVertex, "fnCommonVertex/I");
}
void CommonVertexBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray
  cloneCommonVertex->Clear();
  fnCommonVertex = 0;

  _selectedMuons.clear();
  _selectedElectrons.clear();
  _selectedTaus.clear();

  edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilderHandle);
  const TransientTrackBuilder* builder = trackBuilderHandle.product();

  int indx = 0;
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(_patMuonSrc, muons);
  if (muons.isValid()) {
    for (std::vector<pat::Muon>::const_iterator it  = muons->begin(); 
	                                        it != muons->end(); ++it, indx++) {
      // if not global continue
      if (!it->isGlobalMuon()) continue;
      if (!it->isTrackerMuon()) continue;
      if (it->pt() <= _minPtMuon) continue;
      if (std::abs(it->eta()) > _maxEtaMuon) continue;
      if (it->track()->normalizedChi2() > _maxChi2Muon) continue;

      reco::TrackRef tk = it->innerTrack();
      const reco::HitPattern& hitp = tk->hitPattern(); 
      if (hitp.numberOfValidPixelHits() 
        + hitp.numberOfValidTrackerHits() < _minTrkHitsMuon) continue;

      vhtm::Lepton <pat::Muon> m(&(*it), indx);
      _selectedMuons.push_back(m);
    }
    if (_selectedMuons.size() > 1) 
      std::sort(_selectedMuons.begin(), 
                _selectedMuons.end(), PtComparator<vhtm::Lepton<pat::Muon> >());
  }
  else {
    edm::LogError("CommonVertexBlock") 
      << "Error >> Failed to get pat::Muon Collection for label: " 
      << _patMuonSrc;
  }

  int jndx = 0;
  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(_patElectronSrc, electrons);
  if (electrons.isValid()) {
    for (std::vector<pat::Electron>::const_iterator it  = electrons->begin(); 
	                                            it != electrons->end(); ++it, jndx++) {
      // if not ECAL driven, continue
      if (!it->ecalDrivenSeed()) continue;
      if (!it->gsfTrack().isNonnull()) continue;
      if (it->pt() <= _minPtElectron) continue;
      if (abs(it->eta()) > _maxEtaElectron) continue;

      vhtm::Lepton <pat::Electron> e(&(*it), jndx);
      _selectedElectrons.push_back(e);
    }
    if (_selectedElectrons.size() > 1) 
      std::sort(_selectedElectrons.begin(), 
                _selectedElectrons.end(), PtComparator<vhtm::Lepton<pat::Electron> >());
  }
  else {
    edm::LogError("CommonVertexBlock") 
      << "Error >> Failed to get pat::Electron Collection for label: " 
      << _patElectronSrc;
  }

  int kndx = 0;
  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(_patTauSrc, taus);
  if (taus.isValid()) {
    for (std::vector<pat::Tau>::const_iterator it  = taus->begin(); 
	                                       it != taus->end(); ++it, kndx++) {
      if (it->pt() <= _minPtTau) continue;
      if (it->tauID("decayModeFinding") <= 0.5) continue;
      if (it->tauID("againstMuonTight") <= 0.5) continue;
      if (it->tauID("againstElectronLoose") <= 0.5) continue;
      if (it->tauID("byLooseCombinedIsolationDeltaBetaCorr") <= 0.5) continue;

      vhtm::Lepton <pat::Tau> t(&(*it), kndx);
      _selectedTaus.push_back(t);
    }
    if (_selectedTaus.size() > 1) 
      std::sort(_selectedTaus.begin(), 
                _selectedTaus.end(), PtComparator<vhtm::Lepton<pat::Tau> >());
  }
  else {
    edm::LogError("CommonVertexBlock") 
      << "Error >> Failed to get pat::Tau Collection for label: " 
      << _patTauSrc;
  }
  int nMuon = _selectedMuons.size();
  int nElectron = _selectedElectrons.size();
  int nTau = _selectedTaus.size();

  if (_verbosity > 0)
    std::cout << "nMuon = " << nMuon
              << ", nElectron = " << nElectron
              << ", nTau = " << nTau
              << std::endl;

  // Now form various combinations
  // mu,tau
  vertexB = new ( (*cloneCommonVertex)[0] ) vhtm::CommonVertex();
  if ( nMuon >= 2 && nElectron == 0 && nTau >= 1 ) {
    int iMuon[] = {0,1}; int iElectron[] = {-1}; int iTau[] = {0};
    FitVertex(builder, vertexB, 2, iMuon, 0, iElectron, 1, iTau, "mMt");
  }
  vertexB = new ( (*cloneCommonVertex)[1] ) vhtm::CommonVertex();
  if ( nMuon >= 1 && nElectron == 0 && nTau >= 2 ) {
    int iMuon[] = {0}; int iElectron[1]= {-1}; int iTau[] = {0,1};
    FitVertex(builder, vertexB, 1, iMuon, 0, iElectron, 2, iTau, "mtT");
  }

  // e,tau
  vertexB = new ( (*cloneCommonVertex)[2] ) vhtm::CommonVertex();
  if ( nMuon == 0 && nElectron >= 2 && nTau >= 1 ) {
    int iMuon[] = {-1}; int iElectron[] = {0,1}; int iTau[] = {0};
    FitVertex(builder, vertexB, 0, iMuon, 2, iElectron, 1, iTau, "eEt");
  }
  vertexB = new ( (*cloneCommonVertex)[3] ) vhtm::CommonVertex();
  if ( nMuon == 0 && nElectron >= 1 && nTau >= 2 ) {
    int iMuon[] = {-1}; int iElectron[] = {0}; int iTau[] = {0,1};
    FitVertex(builder, vertexB, 0, iMuon, 1, iElectron, 2, iTau, "etT");
  }

  // e,mu,tau - 3 lepton combinations
  vertexB = new ( (*cloneCommonVertex)[4] ) vhtm::CommonVertex();
  if ( nMuon >= 1 && nElectron >= 1 && nTau >= 1 ) {
    int iMuon[] = {0}; int iElectron[] = {0}; int iTau[] = {0};
    FitVertex(builder, vertexB, 1, iMuon, 1, iElectron, 1, iTau, "met");
  }
  vertexB = new ( (*cloneCommonVertex)[5] ) vhtm::CommonVertex();
  if ( nMuon >= 2 && nElectron >= 1 && nTau >= 1 ) {
    int iMuon[] = {1}; int iElectron[] = {0}; int iTau[] = {0};
    FitVertex(builder, vertexB, 1, iMuon, 1, iElectron, 1, iTau, "Met");
  }
  vertexB = new ( (*cloneCommonVertex)[6] ) vhtm::CommonVertex();
  if ( nMuon >= 1 && nElectron >= 2 && nTau >= 1 ) {
    int iMuon[] = {0}; int iElectron[] = {1}; int iTau[] = {0};
    FitVertex(builder, vertexB, 1, iMuon, 1, iElectron, 1, iTau, "mEt");
  }
  vertexB = new ( (*cloneCommonVertex)[7] ) vhtm::CommonVertex();
  if ( nMuon >= 1 && nElectron >= 1 && nTau >= 2 ) {
    int iMuon[] = {0}; int iElectron[] = {0}; int iTau[] = {1};
    FitVertex(builder, vertexB, 1, iMuon, 1, iElectron, 1, iTau, "meT");
  }

  fnCommonVertex = 8;
}
void CommonVertexBlock::FitVertex(const TransientTrackBuilder* builder, vhtm::CommonVertex* v,
				  int nMuon, int iMuon[], 
                                  int nElectron, int iElectron[], 
                                  int nTau, int iTau[], std::string label) {

  std::vector<reco::TransientTrack> tracks;
  std::vector<int> indices;
  for (int i = 0; i < nMuon; ++i) {
    if (iMuon[i] < 0) continue;
    vhtm::Lepton<pat::Muon> el = _selectedMuons[iMuon[i]];
    const pat::Muon* muon = el.obj();
    reco::TransientTrack transtracks = getTracks(muon, builder);
    if (transtracks.isValid()) {
      tracks.push_back(transtracks);
      indices.push_back(el.indx());
    }
  }
  for (int i = 0; i < nElectron; ++i) {
    if (iElectron[i] < 0) continue;
    vhtm::Lepton<pat::Electron> el = _selectedElectrons[iElectron[i]];
    const pat::Electron* electron = el.obj();
    reco::TransientTrack transtracks = getTracks(electron, builder);
    if (transtracks.isValid()) {
      tracks.push_back(transtracks);
      indices.push_back(el.indx());
    }
  }
  for (int i = 0; i < nTau; ++i) {
    if (iTau[i] < 0) continue;
    vhtm::Lepton<pat::Tau> el = _selectedTaus[iTau[i]];
    const pat::Tau* tau = el.obj();
    reco::TransientTrack transtracks = getTracks(tau, builder);
    if (transtracks.isValid()) {
      tracks.push_back(transtracks);
      indices.push_back(el.indx());
    }
  }

  // Now fit the tracks assuming a common vertex
  if (tracks.size() > 1) {
    KalmanVertexFitter kvf(true);
    TransientVertex vtx = kvf.vertex(tracks);
    v->chi2 = vtx.totalChiSquared();
    v->ndof = vtx.degreesOfFreedom();
    v->label = label;
    for (unsigned int i = 0; i < NEL(v->indices); ++i) {
      v->indices[i] = indices[i];
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CommonVertexBlock);
