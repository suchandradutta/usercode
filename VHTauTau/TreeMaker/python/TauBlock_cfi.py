import FWCore.ParameterSet.Config as cms

from VHTauTau.TreeMaker.AntiElectronIDMVA_cfi import antiElectronMVAIDcfg
tauBlock = cms.EDAnalyzer("TauBlock",
  antiElectronMVAIDcfg,
  verbosity = cms.int32(0),
  patTauSrc = cms.InputTag('selectedPatTaus')
)
