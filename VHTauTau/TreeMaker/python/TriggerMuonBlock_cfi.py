import FWCore.ParameterSet.Config as cms

triggerMuonBlock = cms.EDAnalyzer("TriggerMuonBlock",
  verbosity = cms.int32(0),
  hltInputTag = cms.InputTag('TriggerResults','','HLT'),
  triggerEventTag = cms.InputTag('patTriggerEvent'),
  muonSrcTag = cms.InputTag('selectedPatMuons'),
  hltPathsOfInterest = cms.vstring(
                     'HLT_IsoMu',
                     'HLT_Mu'),
  muonMatchLabel = cms.string( 'matchHLTMuons' ),
  doubleMuonPath = cms.string( 'HLT_Mu13_Mu8_v' )
                                  
)
