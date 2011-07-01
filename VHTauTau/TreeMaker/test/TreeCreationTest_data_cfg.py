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
                      input = cms.untracked.int32(50)
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
process.GlobalTag.globaltag = 'GR_R_42_V14::All'
#from PhysicsTools.PatAlgos.patTemplate_cfg import *
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
     fileName = cms.string('SingleMu.root')
     )                                   
#--------------------------------------------------
# VHTauTau Tree Specific
#--------------------------------------------------
process.load("VHTauTau.TreeMaker.TreeCreator_cfi")
process.load("VHTauTau.TreeMaker.TreeWriter_cfi")
process.load("VHTauTau.TreeMaker.TreeContentConfig_data_cff")

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
                     genJetCollection=cms.InputTag(""),
                     doJetID      = False
                 )

#switchToPFTauHPSpTaNC(process) # For HPS TaNC Taus
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process) # For HPS Taus

# create TC and PF METs
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

process.p = cms.Path(
    process.treeCreator +
    process.PFTau +
    process.patDefaultSequence  +   
    process.treeContentSequence +
    process.treeWriter
    )

from VHTauTau.TreeMaker.SwitchToData import switchToData
switchToData(process)

#--------------------------------------
# List File names here
#---------------------------------------
process.PoolSource.fileNames = [
     '/store/data/Run2011A/SingleMu/AOD/PromptReco-v4/000/165/620/6AD32415-7888-E011-A499-001D09F292D1.root'
]

