import FWCore.ParameterSet.Config as cms


from DQMOffline.PFTau.JetValidation_cfi import JetBenchmark

pfJetValidation1 = JetBenchmark.clone()
pfJetValidation1.InputCollection = cms.InputTag('ak5PFJets')
pfJetValidation1.MatchCollection = cms.InputTag('ak5GenJets')
pfJetValidation1.BenchmarkLabel  = cms.string('PFJet/CompWithGenJet')

pfJetValidation2 = JetBenchmark.clone()
pfJetValidation2.InputCollection = cms.InputTag('ak5PFJets')
pfJetValidation2.MatchCollection = cms.InputTag('ak5CaloJets')
pfJetValidation2.BenchmarkLabel  = cms.string('PFJet/CompWithCaloJet')

pfJetValidationSequence = cms.Sequence( pfJetValidation1 * pfJetValidation2 )
