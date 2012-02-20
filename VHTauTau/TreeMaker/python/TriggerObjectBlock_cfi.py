import FWCore.ParameterSet.Config as cms

triggerObjectBlock = cms.EDAnalyzer("TriggerObjectBlock",
  verbosity = cms.int32(0),
  hltInputTag = cms.InputTag('TriggerResults','','HLT'),
  triggerEventTag = cms.InputTag('patTriggerEvent'),
  hltPathsOfInterest = cms.vstring ("HLT_DoubleMu7_v", "HLT_Mu13_Mu8_v", "HLT_Mu17_Mu8_v", "HLT_Mu8_v",
                                    "HLT_Mu15_v", "HLT_Mu20_v", "HLT_Mu24_v",
                                    "HLT_IsoMu15_v", "HLT_IsoMu20_v", "HLT_IsoMu24_v",
                                    "HLT_Mu17_Ele8_Calo", "HLT_Mu8_Ele17_",
                                    "LooseIsoPFTau", "TightIsoPFTau", "DoubleIsoPFTau"),
  May10ReRecoData = cms.bool(False)                                  
)
