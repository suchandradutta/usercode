import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauTree")
#------------------------
# Message Logger Settings
#------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 200
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
process.GlobalTag.globaltag = 'START42_V13::All'
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
# Event Skim
process.load("VHTauTau.TreeMaker.EventSkimmer_cfi")

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
#-------------------------------------------------------
# PAT 
#------------------------------------------------------
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Add PF jets
from PhysicsTools.PatAlgos.tools.tauTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
jec = [ 'L1FastJet', 'L2Relative', 'L3Absolute' ]
JetCorrectionService = cms.string('ak5PFL1FastL2L3'),
#if not isMC:
#        jec.extend([ 'L2L3Residual' ])
addJetCollection(process, cms.InputTag('ak5PFJets'),
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
     outputModule     = ''
)
addJetCollection(process, cms.InputTag('ak5CaloJets'),
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
     outputModule     = ''
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
           tag    = cms.string('JetCorrectorParametersCollection_Jec11V2_AK5PF'),
           label  = cms.untracked.string('AK5PF')
       ),
       cms.PSet(
           record = cms.string('JetCorrectionsRecord'),
           tag    = cms.string('JetCorrectorParametersCollection_Jec11V2_AK5Calo'),
           label  = cms.untracked.string('AK5Calo')
       )
    ),
    connect = cms.string('sqlite_fip:TauAnalysis/Configuration/data/Jec11V2.db')
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'jec')

##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
##-------------------- Turn-on the FastJet density calculation

process.kt6PFJets.doRhoFastjet = True
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm
process.ak5PFJets.doAreaFastjet = True

# PF based muon and electron isolation
# Muon Isolation
from TauAnalysis.RecoTools.patLeptonPFIsolationSelector_cfi import patMuonPFIsolationSelector
patMuonPFIsolationSelector.chargedHadronIso.ptMin = cms.double(0.5)
patMuonPFIsolationSelector.neutralHadronIso.ptMin = cms.double(0.5)
patMuonPFIsolationSelector.photonIso.ptMin = cms.double(0.5)
patMuonPFIsolationSelector.chargedHadronIso.dRvetoCone = cms.double(0.01)
patMuonPFIsolationSelector.neutralHadronIso.dRvetoCone = cms.double(0.01)
patMuonPFIsolationSelector.photonIso.dRvetoCone = cms.double(0.01)
patMuonPFIsolationSelector.pileUpCorr.chargedToNeutralFactor = cms.double(0.5)
process.patMuonsLoosePFIsoEmbedded04 = cms.EDProducer("PATMuonPFIsolationEmbedder",
    patMuonPFIsolationSelector,
    src = cms.InputTag("selectedPatMuons"),
    userFloatName = cms.string('pfLooseIsoPt04')
)
process.patMuonsLoosePFIsoEmbedded04.pfCandidateSource = cms.InputTag('pfNoPileUp')
process.muonBlock.muonSrc = cms.InputTag("patMuonsLoosePFIsoEmbedded04")

# Electron Isolation
from TauAnalysis.RecoTools.patLeptonPFIsolationSelector_cfi import patElectronPFIsolationSelector
patElectronPFIsolationSelector.chargedHadronIso.ptMin = cms.double(0.5)
patElectronPFIsolationSelector.neutralHadronIso.ptMin = cms.double(0.5)
patElectronPFIsolationSelector.photonIso.ptMin = cms.double(0.5)
patElectronPFIsolationSelector.chargedHadronIso.dRvetoCone = cms.double(0.03)
patElectronPFIsolationSelector.neutralHadronIso.dRvetoCone = cms.double(0.08)
patElectronPFIsolationSelector.photonIso.dRvetoCone = cms.double(0.05)
patElectronPFIsolationSelector.pileUpCorr.chargedToNeutralFactor = cms.double(0.5)
process.patElectronsLoosePFIsoEmbedded04 = cms.EDProducer("PATElectronPFIsolationEmbedder",
    patElectronPFIsolationSelector,
    src = cms.InputTag("selectedPatElectrons"),
    userFloatName = cms.string('pfLooseIsoPt04')
)
process.patElectronsLoosePFIsoEmbedded04.pfCandidateSource = cms.InputTag('pfNoPileUp')
process.electronBlock.electronSrc = cms.InputTag("patElectronsLoosePFIsoEmbedded04")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process) # For HPS Taus
#switchToPFTauHPSpTaNC(process) # For HPS TaNC Taus
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

# remove MC matching
#removeMCMatching(process, ['All'])

process.load("RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff")
#process.patElectronIsolation = cms.Sequence(process.egammaIsolationSequence)

process.patElectrons.isoDeposits = cms.PSet()
process.patElectrons.userIsolation = cms.PSet()
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

process.triggerBlock.hltPathsOfInterest = cms.vstring(
    "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL",
    "HLT_DoubleMu7",
    "HLT_Mu17_Ele8_CaloIdL",
    "HLT_IsoMu12_LooseIsoPFTau10"
)
process.jetBlock.verbosity = cms.int32(1)
process.p = cms.Path(
   process.treeCreator +
   process.kt6PFJets *
   process.ak5PFJets *
   process.PFTau     +
   process.patDefaultSequence +
   process.patMuonsLoosePFIsoEmbedded04 +
   process.patElectronsLoosePFIsoEmbedded04 +
   process.treeContentSequence +
   process.treeWriter
)
#--------------------------------------
# List File names here
#---------------------------------------
process.PoolSource.fileNames = [
     '/store/relval/CMSSW_4_2_6/RelValTTbar/GEN-SIM-RECO/MC_42_V12-v1/0006/AE44BEC1-E8A8-E011-BB78-003048D25B68.root',
     '/store/relval/CMSSW_4_2_6/RelValTTbar/GEN-SIM-RECO/MC_42_V12-v1/0007/203FAE61-3CA9-E011-AA6D-001A928116B8.root'
]
