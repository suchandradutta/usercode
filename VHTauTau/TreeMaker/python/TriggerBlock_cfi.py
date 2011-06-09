import FWCore.ParameterSet.Config as cms

triggerBlock = cms.EDAnalyzer("TriggerBlock",
  verbosity = cms.int32(0),
  l1InputTag  = cms.InputTag('gtDigis'),
  hltInputTag = cms.InputTag('TriggerResults','','HLT'),
  hltPathsOfInterest = cms.vstring(
                                         #Mu
                                         'HLT_Mu12',
                                         'HLT_Mu15',
                                         'HLT_Mu20',
                                         'HLT_Mu24',
                                         'HLT_Mu30',
                                         'HLT_IsoMu12',
                                         'HLT_IsoMu15',
                                         'HLT_IsoMu20',
                                         'HLT_IsoMu24',
                                         'HLT_isoMu30'
                                  )       
)
