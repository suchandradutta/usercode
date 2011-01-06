import FWCore.ParameterSet.Config as cms

MetBenchmark = cms.EDAnalyzer("MetBenchmarkAnalyzer",
    InputCollection = cms.InputTag('pfMet'),
    MatchCollection = cms.InputTag('met'),
    BenchmarkLabel  = cms.string('PFMET/CompWithCaloMET'),
    mode            = cms.int32( 1 ),
    VariableEtBins  = cms.vdouble(0.,1.,2.,5.,10.,20.,50.,100.,200.,400.,1000.),
    CreateMETSpecificHistos = cms.bool(True),
    ptMin = cms.double( 0.0 ),
    ptMax = cms.double( 999999 ),
    etaMin = cms.double(-10),
    etaMax = cms.double(10),
    phiMin = cms.double(-3.14),
    phiMax = cms.double(3.14)
)
