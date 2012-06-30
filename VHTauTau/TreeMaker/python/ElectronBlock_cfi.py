import FWCore.ParameterSet.Config as cms

from VHTauTau.TreeMaker.ElectronIDMVA_cfi import electronMVAIDNOIPcfg
from VHTauTau.TreeMaker.ElectronMVA_cfi import electronIdMVAcfg,electronIsoMVAcfg

electronBlock = cms.EDAnalyzer("ElectronBlock",
  electronMVAIDNOIPcfg,
  electronIdMVAcfg,
  electronIsoMVAcfg,
  verbosity       = cms.int32(0),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  trackSrc        = cms.InputTag('generalTracks'),
  dcsSrc          = cms.InputTag('scalersRawToDigi'),
  vertexSrc       = cms.InputTag('offlinePrimaryVerticesWithBS'),
  convSrc         = cms.InputTag('allConversions'),
  electronSrc     = cms.InputTag('selectedPatElectrons'),
  pfElectronSrc   = cms.InputTag('particleFlow'),
  ecalEBSrc       = cms.InputTag('reducedEcalRecHitsEB'),
  ecalEESrc       = cms.InputTag('reducedEcalRecHitsEE'),
  rhoSrc          = cms.InputTag('kt6PFJets','rho'),
  pfSrc           = cms.InputTag('particleFlow')
)
