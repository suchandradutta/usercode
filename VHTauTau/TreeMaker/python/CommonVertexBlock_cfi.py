import FWCore.ParameterSet.Config as cms

commonVertexBlock = cms.EDAnalyzer("CommonVertexBlock",
  verbosity = cms.int32(0),
  patMuonSrc = cms.InputTag('cleanPatMuons'),
  patElectronSrc = cms.InputTag('cleanPatElectrons'),
  patTauSrc = cms.InputTag('cleanPatTaus'),
  minPtMuon = cms.double(8.0),
  maxEtaMuon = cms.double(2.3),
  maxChi2Muon = cms.double(15.),
  minTrkHitsMuon = cms.double(8.),
  minPtElectron = cms.double(8.0),
  maxEtaElectron = cms.double(2.3),
  minPtTau = cms.double(10.0)
)
