import FWCore.ParameterSet.Config as cms

caloJetBlock       = cms.EDAnalyzer("CaloJetBlock",
    verbosity      = cms.int32(0),
    caloJetSrc     = cms.InputTag('cleanPatJets'),
    electronPt     = cms.double(30.),
    electronIso    = cms.double(0.1), 
    muonPt         = cms.double(10.),
    muonIso        = cms.double(0.05),
    jecUncertainty = cms.string('CondFormats/JetMETObjects/data/Spring10_Uncertainty_AK5Calo.txt'),
    applyResidualJEC = cms.bool(False),
    residualJEC = cms.string('CondFormats/JetMETObjects/data/Spring10_L2L3Residual_AK5Calo.txt')
)
