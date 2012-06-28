import FWCore.ParameterSet.Config as cms

SVfitBlock = cms.EDAnalyzer("SVfitBlock",
  verbosity    = cms.int32(0),
  diTauPairs   = cms.InputTag('DiTauPairsWHInt'),
  genParticles = cms.InputTag('genParticles')
)
