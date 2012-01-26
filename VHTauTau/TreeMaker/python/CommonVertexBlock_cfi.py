import FWCore.ParameterSet.Config as cms

commonVertexBlock = cms.EDAnalyzer("CommonVertexBlock",
  verbosity = cms.int32(0),
  patMuonSrc = cms.InputTag('cleanPatMuons'),
  patElectronSrc = cms.InputTag('cleanPatElectrons'),
  patTauSrc = cms.InputTag('cleanPatTaus'),
  minPtMuon = cms.double(5.0),
  minPtElectron = cms.double(5.0),
  minPtTau = cms.double(10.0)
)
