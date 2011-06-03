import FWCore.ParameterSet.Config as cms

electronBlock      = cms.EDAnalyzer("ElectronBlock",
    verbosity      = cms.int32(0),
    trackSrc       = cms.InputTag('generalTracks'),
    dcsSrc         = cms.InputTag('scalersRawToDigi'),
    vertexSrc      = cms.InputTag('offlinePrimaryVertices'),
    electronSrc    = cms.InputTag('cleanPatElectrons'),
    pfElectronSrc  = cms.InputTag('particleFlow')                        
)
