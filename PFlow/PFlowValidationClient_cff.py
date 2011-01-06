import FWCore.ParameterSet.Config as cms

from DQMOffline.PFTau.BenchmarkClient_cfi import benchmarkClient

pfJetClient = benchmarkClient.clone()
pfJetClient.FolderNames = cms.vstring("PFJet/CompWithGenJet","PFJet/CompWithCaloJet")
pfJetClient.HistogramNames = cms.vstring( "delta_et_VS_et_")

pfMETClient = benchmarkClient.clone()
pfMETClient.FolderNames = cms.vstring("PFMET/CompWithCaloMET")
pfMETClient.HistogramNames = cms.vstring( "delta_et_VS_et_")
