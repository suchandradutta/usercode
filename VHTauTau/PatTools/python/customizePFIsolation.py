import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFMuonIso, setupPFElectronIso
from CommonTools.ParticleFlow.pfParticleSelection_cff import pfParticleSelectionSequence

def configurePFIsolation(process):
  process.muIsoSequence       = setupPFMuonIso(process,'muons')
  process.electronIsoSequence = setupPFElectronIso(process,'gsfElectrons')
  process.pfParticleSelectionSequence = pfParticleSelectionSequence

  process.patMuons.isoDeposits = cms.PSet(
    pfAllParticles   = cms.InputTag("muPFIsoDepositPUPFIso"),      # all PU   CH+MU+E
    pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedPFIso"), # all noPU CH
    pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralPFIso"), # all NH
    pfPhotons        = cms.InputTag("muPFIsoDepositGammaPFIso"),   # all PH
    user = cms.VInputTag(
      cms.InputTag("muPFIsoDepositChargedAllPFIso")                # all noPU CH+MU+E
    )
  )
  process.patMuons.isolationValues = cms.PSet(
    pfAllParticles   = cms.InputTag("muPFIsoValuePU04PFIso"),
    pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04PFIso"),
    pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04PFIso"),
    pfPhotons        = cms.InputTag("muPFIsoValueGamma04PFIso"),
    user = cms.VInputTag(
      cms.InputTag("muPFIsoValueChargedAll04PFIso")
    )
  )
  process.patElectrons.isoDeposits = cms.PSet(
    pfAllParticles   = cms.InputTag("elPFIsoDepositPUPFIso"),      # all PU   CH+MU+E
    pfChargedHadrons = cms.InputTag("elPFIsoDepositChargedPFIso"), # all noPU CH
    pfNeutralHadrons = cms.InputTag("elPFIsoDepositNeutralPFIso"), # all NH
    pfPhotons        = cms.InputTag("elPFIsoDepositGammaPFIso"),   # all PH
    user = cms.VInputTag(
      cms.InputTag("elPFIsoDepositChargedAllPFIso"),                 # all noPU CH+MU+E
    )
  )
  process.patElectrons.isolationValues = cms.PSet(
    pfAllParticles   = cms.InputTag("elPFIsoValuePU04PFIdPFIso"),
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04PFIdPFIso"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04PFIdPFIso"),
    pfPhotons        = cms.InputTag("elPFIsoValueGamma04PFIdPFIso"),
    user = cms.VInputTag(
      cms.InputTag("elPFIsoValueChargedAll04PFIdPFIso"),
      cms.InputTag("elPFIsoValueChargedAll04NoPFIdPFIso"),
      cms.InputTag("elPFIsoValuePU04NoPFIdPFIso"),
      cms.InputTag("elPFIsoValueCharged04NoPFIdPFIso"),
      cms.InputTag("elPFIsoValueGamma04NoPFIdPFIso"),
      cms.InputTag("elPFIsoValueNeutral04NoPFIdPFIso")
    )
  )
  process.pfIsolationSequence = cms.Sequence(
    process.pfParticleSelectionSequence +
    process.process.muIsoSequence +
    process.electronIsoSequence
  )
  return process.pfIsolationSequence
