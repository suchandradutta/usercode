import FWCore.ParameterSet.Config as cms

superClusterBlock = cms.EDAnalyzer("SuperClusterBlock",
  verbosity         = cms.int32(0),
  ebInputTag        = cms.InputTag('correctedHybridSuperClusters'),
  eeInputTag        = cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'),
  ecalEBInputTag    = cms.InputTag('reducedEcalRecHitsEB'),
  ecalEEInputTag    = cms.InputTag('reducedEcalRecHitsEE'),
  tracksInputTag    = cms.InputTag('generalTracks'),
  electronsInputTag = cms.InputTag('selectedPatElectrons')
)
