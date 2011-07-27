import FWCore.ParameterSet.Config as cms

photonBlock = cms.EDAnalyzer("PhotonBlock",
  verbosity = cms.int32(0),
  treeName  = cms.string('vhtree'),
  photonSrc = cms.InputTag('cleanPatPhotons') 
)
