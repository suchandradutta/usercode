import FWCore.ParameterSet.Config as cms

JetBenchmark = cms.EDAnalyzer("JetBenchmarkAnalyzer",
    InputCollection = cms.InputTag('ak5PFJets'),
    MatchCollection = cms.InputTag('ak5CaloJets'),
    BenchmarkLabel  = cms.string('ParticleFlow/PFVsCalo'),
    deltaRMax = cms.double(0.1),
    matchCharge = cms.bool(False),
    mode = cms.int32( 1 ),
    CreatePFractionHistos = cms.bool(False),
    ptMin = cms.double( 0.0 ),
    ptMax = cms.double( 999999 ),
    etaMin = cms.double(-10),
    etaMax = cms.double(10),
    phiMin = cms.double(-3.14),
    phiMax = cms.double(3.14)
)
