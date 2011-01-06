import FWCore.ParameterSet.Config as cms


from DQMOffline.PFTau.MetValidation_cfi import MetBenchmark

pfMetValidation1 = MetBenchmark.clone()
pfMetValidation1.InputCollection = cms.InputTag('pfMet')
pfMetValidation1.MatchCollection = cms.InputTag('genMet')
pfMetValidation1.BenchmarkLabel  = cms.string('PFMET/CompWithGenMET')

pfMetValidation2 = MetBenchmark.clone()
pfMetValidation2.InputCollection = cms.InputTag('pfMet')
pfMetValidation2.MatchCollection = cms.InputTag('met')
pfMetValidation2.BenchmarkLabel  = cms.string('PFMET/CompWithCaloMET')

pfMetValidationSequence = cms.Sequence( pfMetValidation1 * pfMetValidation2 )
