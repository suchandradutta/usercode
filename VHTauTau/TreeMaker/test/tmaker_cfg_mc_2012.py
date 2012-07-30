import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauTree")
#------------------------
# Message Logger Settings
#------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring()
)
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(100)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

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
process.GlobalTag.globaltag = 'START52_V11::All'
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('tree.root')
)
#--------------------------------------------------
# VHTauTau Tree Specific
#--------------------------------------------------
process.load("VHTauTau.TreeMaker.TreeCreator_cfi")
process.load("VHTauTau.TreeMaker.TreeWriter_cfi")
process.load("VHTauTau.TreeMaker.TreeContentConfig_data_cff")
#-------------------------------------------------------
# PAT 
#------------------------------------------------------
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")

import PhysicsTools.PatAlgos.tools.tauTools as tauTools
import PhysicsTools.PatAlgos.tools.jetTools as jetTools
import PhysicsTools.PatAlgos.tools.metTools as metTools
import PhysicsTools.PatAlgos.tools.pfTools as pfTools

tauTools.switchToPFTauHPS(process) # For HPS Taus

metTools.addTcMET(process, 'TC')
metTools.addPfMET(process, 'PF')

## --
## Switch on PAT trigger
## --
import PhysicsTools.PatAlgos.tools.trigTools as trigTools
trigTools.switchOnTrigger( process, outputModule='' ) # This is optional and can be omitted.

# load the PU JetID sequence
process.load("CMGTools.External.pujetidsequence_cff")
process.puJetId.jets = cms.InputTag("selectedPatJetsAK5PF")
process.puJetMva.jets = cms.InputTag("selectedPatJetsAK5PF")

jec = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
jetTools.addJetCollection(process, cms.InputTag('ak5PFJets'),
  'AK5', 'PF',
   doJTA            = True,
   doBTagging       = True,
   jetCorrLabel     = ('AK5PF', cms.vstring(jec)),
   doType1MET       = False,
   doL1Cleaning     = True,
   doL1Counters     = False,
   genJetCollection = cms.InputTag("ak5GenJets"),
   doJetID          = True,
   jetIdLabel       = "ak5",
   outputModules    = ['']
)
jetTools.addJetCollection(process, cms.InputTag('ak5CaloJets'),
   'AK5', 'Calo',
   doJTA            = True,
   doBTagging       = True,
   jetCorrLabel     = ('AK5Calo', cms.vstring(jec)),
   doType1MET       = True,
   doL1Cleaning     = True,
   doL1Counters     = False,
   genJetCollection = cms.InputTag("ak5GenJets"),
   doJetID          = True,
   jetIdLabel       = "ak5",
   outputModules    = ['']
)
#--------------------------------------------------------------------------------
#
# configure Jet Energy Corrections
#
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.jec = cms.ESSource("PoolDBESSource",
  DBParameters = cms.PSet(
    messageLevel = cms.untracked.int32(0)
  ),
  timetype = cms.string('runnumber'),
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK5PF'),
      label  = cms.untracked.string('AK5PF')
    ),
    cms.PSet(
      record = cms.string('JetCorrectionsRecord'),
      tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK5Calo'),
      label  = cms.untracked.string('AK5Calo')
    )
  ),
  connect = cms.string('sqlite_fip:VHTauTau/TreeMaker/data/Summer12_V7_MC.db')
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')
#-------------------------------------------------------------------------------------------------------------------------
##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation
process.kt6PFJets.doRhoFastjet = True
process.kt6PFJets.Rho_EtaMax = cms.double(4.4)
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm 
process.ak5PFJets.doAreaFastjet = True
process.ak5PFJets.Rho_EtaMax = cms.double(4.4)
#process.patJetCorrFactors.useRho = True

## re-run kt4PFJets within lepton acceptance to compute rho
process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJetsCentral = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsCentral.Rho_EtaMax = cms.double(2.5)

process.fjSequence = cms.Sequence(process.kt6PFJets+process.ak5PFJets+process.kt6PFJetsCentral)

process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff")
#process.patElectronIsolation = cms.Sequence(process.egammaIsolationSequence)

#process.patElectrons.isoDeposits = cms.PSet()
#process.patElectrons.userIsolation = cms.PSet()
process.patElectrons.addElectronID = cms.bool(True)
process.patElectrons.electronIDSources = cms.PSet(
    simpleEleId95relIso = cms.InputTag("simpleEleId95relIso"),
    simpleEleId90relIso = cms.InputTag("simpleEleId90relIso"),
    simpleEleId85relIso = cms.InputTag("simpleEleId85relIso"),
    simpleEleId80relIso = cms.InputTag("simpleEleId80relIso"),
    simpleEleId70relIso = cms.InputTag("simpleEleId70relIso"),
    simpleEleId60relIso = cms.InputTag("simpleEleId60relIso"),
    simpleEleId95cIso = cms.InputTag("simpleEleId95cIso"),
    simpleEleId90cIso = cms.InputTag("simpleEleId90cIso"),
    simpleEleId85cIso = cms.InputTag("simpleEleId85cIso"),
    simpleEleId80cIso = cms.InputTag("simpleEleId80cIso"),
    simpleEleId70cIso = cms.InputTag("simpleEleId70cIso"),
    simpleEleId60cIso = cms.InputTag("simpleEleId60cIso")
)
##
process.patElectrons.addGenMatch = cms.bool(False)
process.patElectrons.embedGenMatch = cms.bool(False)
#process.patElectrons.usePV = cms.bool(False)
##
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
# you have to tell the ID that it is data
process.simpleEleId95relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId90relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId85relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId80relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId70relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId60relIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId95cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId90cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId85cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId80cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId70cIso.dataMagneticFieldSetUp = cms.bool(True)
process.simpleEleId60cIso.dataMagneticFieldSetUp = cms.bool(True)
#
process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.makePatElectrons = cms.Sequence(process.patElectronIDs*process.patElectrons) 

# MVA MET
process.load("RecoMET/METProducers/mvaPFMET_cff")
process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring("ak5PFL1FastL2L3Residual")
#process.pfMEtMVA.srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatTaus")
process.pfMEtMVA.srcLeptons = cms.VInputTag('cleanPatElectrons', 'cleanPatMuons', 'cleanPatTaus')

process.patPFMetByMVA = process.patMETs.clone(
  metSource = cms.InputTag('pfMEtMVA'),
  addMuonCorrections = cms.bool(False),
  addGenMET = cms.bool(False),
  genMETSource = cms.InputTag('genMetTrue')
)

# SVfit
from VHTauTau.TreeMaker.customizeSVfit import configureSVfit
process.SVND = configureSVfit(process)

# pfIsolation
from VHTauTau.PatTools.customizePFIsolation import configurePFIsolation
process.pfIsolationSeq = configurePFIsolation(process)

# Disable tau IsoDeposits
process.patDefaultSequence.remove(process.patPFTauIsolation)
process.patTaus.isoDeposits = cms.PSet()
process.patTaus.userIsolation = cms.PSet()

process.p = cms.Path(
  process.fjSequence +
#  process.PFTau +
  process.pfIsolationSeq +  
  process.patDefaultSequence +
  process.puJetIdSqeuence + 
  process.pfMEtMVAsequence + 
  process.patPFMetByMVA +
  process.SVND + 
  process.treeCreator +
  process.treeContentSequence +
  process.treeWriter
)

# List File names here
#---------------------------------------
process.PoolSource.fileNames = [
  '/store/mc/Summer12/WH_ZH_TTH_HToTauTau_M-125_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v2/0000/0A439819-32A5-E111-9129-00215E22190E.root'
]