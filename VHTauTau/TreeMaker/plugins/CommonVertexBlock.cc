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

namespace {
  template <class T>
  class PtComparator {
  public:
    bool operator()(const T* a, const T* b) const {
      return a->pt() > b->pt();
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
  _minPtElectron(iConfig.getParameter<double>("minPtElectron")),
  _minPtTau(iConfig.getParameter<double>("minPtTau"))
{}
void CommonVertexBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneCommonVertex = new TClonesArray("vhtm::CommonVertex");
  tree->Branch("CommonVertex", &cloneCommonVertex, 32000, 2);
}
void CommonVertexBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray
  cloneCommonVertex->Clear();

  _selectedMuons.clear();
  _selectedElectrons.clear();
  _selectedTaus.clear();

  edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilderHandle);
  const TransientTrackBuilder* builder = trackBuilderHandle.product();

  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(_patMuonSrc, muons);
  if (muons.isValid()) {
    for (std::vector<pat::Muon>::const_iterator it  = muons->begin(); 
                                                it != muons->end(); ++it) {
      // if not global continue
      if (!it->isGlobalMuon()) continue;
      if (it->pt() <= _minPtMuon) continue;

      _selectedMuons.push_back(&(*it));
    }
    if (_selectedMuons.size() > 1) 
      std::sort(_selectedMuons.begin(), 
                _selectedMuons.end(), PtComparator<pat::Muon>());
  }
  else {
    edm::LogError("CommonVertexBlock") 
      << "Error! Failed to get PAT Muon Collection for tag = " 
      << _patMuonSrc;
  }

  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(_patElectronSrc, electrons);
  if (electrons.isValid()) {
    for (std::vector<pat::Electron>::const_iterator it  = electrons->begin(); 
                                                    it != electrons->end(); ++it) {
      // if not ECAL driven, continue
      if (!it->ecalDrivenSeed()) continue;
      if (it->pt() <= _minPtElectron) continue;

      _selectedElectrons.push_back(&(*it));
    }
    if (_selectedElectrons.size() > 1) 
      std::sort(_selectedElectrons.begin(), 
                _selectedElectrons.end(), PtComparator<pat::Electron>());
  }
  else {
    edm::LogError("CommonVertexBlock") 
      << "Error! Failed to get PAT Electron Collection for tag = " 
      << _patElectronSrc;
  }

  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(_patTauSrc, taus);
  if (taus.isValid()) {
    for (std::vector<pat::Tau>::const_iterator it  = taus->begin(); 
                                               it != taus->end(); ++it) {
      if (it->pt() <= _minPtTau) continue;
      if (it->tauID("decayModeFinding") <= 0.5) continue;

      _selectedTaus.push_back(&(*it));
    }
    if (_selectedTaus.size() > 1) 
      std::sort(_selectedTaus.begin(), 
                _selectedTaus.end(), PtComparator<pat::Tau>());
  }
  else {
    edm::LogError("CommonVertexBlock") 
      << "Error! Failed to get PAT Tau Collection for tag = " 
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
  vertexB = new ( (*cloneCommonVertex)[0] ) vhtm::CommonVertex();
  if ( nMuon >= 2 && nElectron == 0 && nTau >= 1 ) 
    FitVertex(builder, vertexB, 2, 0, 1, "mmt");

  vertexB = new ( (*cloneCommonVertex)[1] ) vhtm::CommonVertex();
  if ( nMuon >= 1 && nElectron == 0 && nTau >= 2 ) 
    FitVertex(builder, vertexB, 1, 0, 2, "mtt");

  vertexB = new ( (*cloneCommonVertex)[2] ) vhtm::CommonVertex();
  if ( nMuon == 0 && nElectron >= 2 && nTau >= 1 ) 
    FitVertex(builder, vertexB, 0, 2, 1, "eet");

  vertexB = new ( (*cloneCommonVertex)[3] ) vhtm::CommonVertex();
  if ( nMuon == 0 && nElectron >= 1 && nTau >= 2 ) 
    FitVertex(builder, vertexB, 0, 1, 2, "ett");

  vertexB = new ( (*cloneCommonVertex)[4] ) vhtm::CommonVertex();
  if ( nMuon >= 1 && nElectron >= 1 && nTau >= 1 ) 
    FitVertex(builder, vertexB, 1, 1, 1, "emt");
}
void CommonVertexBlock::FitVertex(const TransientTrackBuilder* builder, vhtm::CommonVertex* v, 
                             int nMuon, int nElectron, int nTau, std::string label) {
  std::vector<reco::TransientTrack> tracks;
  for (int i = 0; i < nMuon; ++i) {
    const pat::Muon* muon = _selectedMuons[i];
    reco::TransientTrack transtracks = getTracks(muon, builder);
    if (transtracks.isValid()) tracks.push_back(transtracks); 
  }
  for (int i = 0; i < nElectron; ++i) {
    const pat::Electron* electron = _selectedElectrons[i];
    reco::TransientTrack transtracks = getTracks(electron, builder);
    if (transtracks.isValid()) tracks.push_back(transtracks); 
  }
  for (int i = 0; i < nTau; ++i) {
    const pat::Tau* tau = _selectedTaus[i];
    reco::TransientTrack transtracks = getTracks(tau, builder);
    if (transtracks.isValid()) tracks.push_back(transtracks); 
  }

  // Now fit the tracks assuming a common vertex
  if (tracks.size() > 1) {
    KalmanVertexFitter kvf(true);
    TransientVertex vtx = kvf.vertex(tracks);
    v->chi2 = vtx.totalChiSquared();
    v->ndof = vtx.degreesOfFreedom();
    v->label = label;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CommonVertexBlock);
