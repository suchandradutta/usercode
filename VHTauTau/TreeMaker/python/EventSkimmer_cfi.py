import FWCore.ParameterSet.Config as cms

eventSkimmer = cms.EDFilter("EventSkimmer",
    verbosity           = cms.int32(0),
    muonSrc             = cms.InputTag('muons'),
    electronSrc         = cms.InputTag('pixelMatchGsfElectrons'),
    doMuonSelection     = cms.bool(True),
    doElectronSelection = cms.bool(False),
    muonPtCut           = cms.double(8.0),
    electronPtCut       = cms.double(8.0),
    createHisto         = cms.bool(True)                            
)
