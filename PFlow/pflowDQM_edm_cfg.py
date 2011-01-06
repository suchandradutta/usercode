# test file for PFCandidate validation
# performs a matching with the genParticles collection. 
# creates a root file with histograms filled with PFCandidate data,
# present in the Candidate, and in the PFCandidate classes, for matched
# PFCandidates. Matching histograms (delta pt etc) are also available. 

import FWCore.ParameterSet.Config as cms

process = cms.Process("PFlowDQM")
process.load("DQMServices.Core.DQM_cfg")

fa = 'RelValQCD'
fb = 'FlatPt_15_3000_Fast'
fc = 'ParticleFlow'

process.source = cms.Source("PoolSource",
         fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_9_5/RelValQCD_FlatPt_15_3000/GEN-SIM-DIGI-RECO/MC_39Y_V6_FastSim-v1/0008/266E2E4B-89FA-DF11-ABA2-0018F3D09710.root'
    )
)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MC_39Y_V6::All'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load("Validation.RecoParticleFlow.JetValidation_cff")
process.load("Validation.RecoParticleFlow.MetValidation_cff")

process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.EDM = cms.OutputModule("PoolOutputModule",
               outputCommands = cms.untracked.vstring('drop *',"keep *_MEtoEDMConverter_*_PFlowDQM"),
               fileName = cms.untracked.string('MEtoEDM_PFlow.root')
               )

process.p =cms.Path(
    process.pfJetValidationSequence +
    process.pfMetValidationSequence +
    process.MEtoEDMConverter
    )


process.outpath = cms.EndPath(process.EDM)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
