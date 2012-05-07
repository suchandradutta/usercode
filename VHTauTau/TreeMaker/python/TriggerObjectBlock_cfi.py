import FWCore.ParameterSet.Config as cms

triggerObjectBlock = cms.EDAnalyzer("TriggerObjectBlock",
  verbosity = cms.int32(0),
  hltInputTag = cms.InputTag('TriggerResults','','HLT'),
  triggerEventTag = cms.InputTag('patTriggerEvent'),
  hltPathsOfInterest = cms.vstring ("HLT_DoubleMu",
                                    "HLT_Mu13_Mu8_v",
                                    "HLT_Mu17_Mu8_v",
                                    "HLT_Mu8_v",
                                    "HLT_Mu1",
                                    "HLT_Mu2",
                                    "HLT_Mu3",
                                    "HLT_Mu4",
                                    "HLT_IsoMu1",
                                    "HLT_IsoMu2",
                                    "HLT_IsoMu3",
                                    "HLT_IsoMu4",
                                    "HLT_Mu17_Ele8_Calo",
                                    "HLT_Mu8_Ele17_",
                                    "IsoPFTau"),
  May10ReRecoData = cms.bool(False)
)
