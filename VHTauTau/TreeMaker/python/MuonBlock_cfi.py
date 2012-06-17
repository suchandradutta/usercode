import FWCore.ParameterSet.Config as cms

from VHTauTau.TreeMaker.MuonMVA_cfi import muonIdMVAcfg,muonIsoMVAcfg,muonIsoRingsRadMVAcfg
muonBlock = cms.EDAnalyzer("MuonBlock",
  muonIdMVAcfg,
  muonIsoMVAcfg,
  muonIsoRingsRadMVAcfg,
  verbosity       = cms.int32(0),
  muonSrc         = cms.InputTag('selectedPatMuons'),
  pfMuonSrc       = cms.InputTag('particleFlow'),
  vertexSrc       = cms.InputTag('offlinePrimaryVerticesWithBS'),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  beamSpotCorr    = cms.bool(True),
  muonID          = cms.string('GlobalMuonPromptTight'),
  rhoSrc          = cms.InputTag('kt6PFJets','rho'),
  pfSrc           = cms.InputTag('particleFlow')
)
