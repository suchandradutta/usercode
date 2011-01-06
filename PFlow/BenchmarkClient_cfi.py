import FWCore.ParameterSet.Config as cms

benchmarkClient = cms.EDAnalyzer("BenchmarkClient",
    FolderNames = cms.vstring("PFJet/CompWithGenJet","PFJet/CompWithCaloJet"),
    HistogramNames = cms.vstring( "delta_et_VS_et_")                                 
)
