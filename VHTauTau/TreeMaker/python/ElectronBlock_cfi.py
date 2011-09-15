import FWCore.ParameterSet.Config as cms

electronBlock = cms.EDAnalyzer("ElectronBlock",
    verbosity     = cms.int32(0),
    trackSrc      = cms.InputTag('generalTracks'),
    dcsSrc        = cms.InputTag('scalersRawToDigi'),
    vertexSrc     = cms.InputTag('offlinePrimaryVertices'),
    electronSrc   = cms.InputTag('selectedPatElectrons'),
    pfElectronSrc = cms.InputTag('particleFlow'),
    ecalEBInputTag = cms.InputTag('reducedEcalRecHitsEB'),
    ecalEEInputTag = cms.InputTag('reducedEcalRecHitsEE')
)
