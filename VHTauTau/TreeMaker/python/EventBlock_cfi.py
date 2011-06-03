import FWCore.ParameterSet.Config as cms

eventBlock         = cms.EDAnalyzer("EventBlock",
    verbosity      = cms.int32(0)
)
