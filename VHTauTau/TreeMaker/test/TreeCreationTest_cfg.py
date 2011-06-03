import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauTree")
#------------------------
# Message Logger Settings
#------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------
process.source = cms.Source("PoolSource",
                   fileNames = cms.untracked.vstring()
                 )
process.maxEvents = cms.untracked.PSet(
                      input = cms.untracked.int32(-1)
                    )
#-----------------------------
# Geometry
#-----------------------------
process.load("Configuration.StandardSequences.Geometry_cff")
#-----------------------------
# Magnetic Field
#-----------------------------
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#-------------
# Global Tag
#-------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START42_V12::All'
#from PhysicsTools.PatAlgos.patTemplate_cfg import *
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
     fileName = cms.string('MC_RelValTTbar.root')
)
#--------------------------------------------------
# VHTauTau Tree Specific
#--------------------------------------------------
process.load("VHTauTau.TreeMaker.TreeCreator_cfi")
process.load("VHTauTau.TreeMaker.TreeWriter_cfi")
process.load("VHTauTau.TreeMaker.TreeContentConfig_cff")

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
#-------------------------------------------------------
# PAT 
#------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Add PF jets
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5PFJets'),
                     'AK5', 'PF',
                     doJTA        = True,
                     doBTagging   = True,
                     jetCorrLabel = ('AK5PF', cms.vstring(['L2Relative', 'L3Absolute', 'L5Flavor', 'L7Parton'])),
                     doType1MET   = False,
                     doL1Cleaning = False,
                     doL1Counters = False,
                     genJetCollection=cms.InputTag("ak5GenJets"),
                     doJetID      = False
                 )

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process) # For HPS Taus
#switchToPFTauHPSpTaNC(process) # For HPS TaNC Taus
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

# remove MC matching
removeMCMatching(process, ['All'])


process.p = cms.Path(
    process.treeCreator +
    process.patDefaultSequence +
    process.PFTau     +
    process.treeContentSequence +
    process.treeWriter
    )

#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('patTuple.root'),
#                               # save only events passing the full path
#                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
#                               # save PAT Layer 1 output; you need a '*' to
#                               # unpack the list of commands 'patEventContent'
#                               outputCommands = cms.untracked.vstring('drop *', *patEventContent )
#                               )

#process.outpath = cms.EndPath(process.out)
#--------------------------------------
# List File names here
#---------------------------------------
process.PoolSource.fileNames = [
#     'rfio:/castor/cern.ch/user/l/lusito/WH/TestA/WH115DIGI2RAWExt.root'
     '/store/relval/CMSSW_4_2_3/RelValTTbar//GEN-SIM-RECO/START42_V12-v2/0068/30222246-647C-E011-A6A6-00304867C1BC.root'    
]

