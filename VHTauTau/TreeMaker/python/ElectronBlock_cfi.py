import FWCore.ParameterSet.Config as cms

from VHTauTau.TreeMaker.ElectronIDMVA_cfi import electronMVAIDNOIPcfg

electronBlock = cms.EDAnalyzer("ElectronBlock",
    electronMVAIDNOIPcfg,
    verbosity     = cms.int32(0),
    offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
    trackSrc      = cms.InputTag('generalTracks'),
    dcsSrc        = cms.InputTag('scalersRawToDigi'),
    vertexSrc     = cms.InputTag('offlinePrimaryVerticesWithBS'),
    convSrc       = cms.InputTag('allConversions'),
    electronSrc   = cms.InputTag('selectedPatElectrons'),
    pfElectronSrc = cms.InputTag('particleFlow'),
    ecalEBInputTag = cms.InputTag('reducedEcalRecHitsEB'),
    ecalEEInputTag = cms.InputTag('reducedEcalRecHitsEE')
)
