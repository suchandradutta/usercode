import FWCore.ParameterSet.Config as cms

import copy
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.tauSelector_cfi import *
from CommonTools.ParticleFlow.TopProjectors.pfNoJet_cfi import *
from RecoMET.METProducers.METSigParams_cfi import *

def configureSVfit(process):
  process.load("TauAnalysis.CandidateTools.nSVfitAlgorithmDiTau_cfi")
  process.load("TauAnalysis.CandidateTools.nSVfitAlgorithmWH_cfi")
  #process.load("PhysicsTools.PatAlgos.patSequences_cff")

  process.selectedPatTausForSVfit = process.selectedPatTaus.clone(
    cut = cms.string('pt > 20 & abs(eta) < 2.3 & tauID("decayModeFinding") > 0.5 & tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 & tauID("againstMuonLoose") > 0.5 & tauID("againstElectronLoose") > 0.5')
  )
  process.selectedPatMuonsForSVfit = process.selectedPatMuons.clone(
    cut = cms.string('pt > 15 & abs(eta) < 2.1 & isGlobalMuon()')
  )

  #-------select ak5PFJets not overlapping with leptons--------
  process.ak5PFJetsNotOverlappingWithLeptons = cms.EDFilter("PFJetAntiOverlapSelector",
    src = cms.InputTag('ak5PFJets'),
    srcNotToBeFiltered = cms.VInputTag(
      'patElectrons',
      'patMuons',
      'patTaus'
    ),
    dRmin = cms.double(0.5),
    invert = cms.bool(False),
    filter = cms.bool(False)
  )

  #--------select PFCandidates ("unclustered energy") not within jets for Type 2 MET correction----------
  process.pfCandsNotInJet = pfNoJet.clone(
    topCollection = 'ak5PFJets',
    bottomCollection = 'particleFlow'
  )

  #--------produce PFMET significance cov. matrix----------
  process.pfMEtSignCovMatrix = cms.EDProducer("PFMEtSignCovMatrixProducer",
    METSignificance_params,
    src = cms.VInputTag(
      'patElectrons',
      'patMuons',
      'patTaus',
      'ak5PFJetsNotOverlappingWithLeptons',
      'pfCandsNotInJet'
    ),
    addJERcorr = cms.PSet(
      inputFileName = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
      lutName = cms.string('pfJetResolutionMCtoDataCorrLUT')
    )
  )

  #------For WH tau+tau SVfit in Fit mode---------
  process.DiTauPairsWHFit = copy.deepcopy(process.nSVfitProducerByLikelihoodMaximizationWH)
  process.DiTauPairsWHFit.config.event.likelihoodFunctions[0].srcMEtCovMatrix = cms.InputTag('pfMEtSignCovMatrix')
  process.DiTauPairsWHFit.config.event.resonances.W.builder.daughters.chargedLepton.src = cms.InputTag("selectedPatMuonsForSVfit")
  process.DiTauPairsWHFit.config.event.resonances.W.daughters.chargedLepton.src = cms.InputTag("selectedPatMuonsForSVfit")
  process.DiTauPairsWHFit.config.event.resonances.Higgs.daughters.leg1.src = cms.InputTag("selectedPatTausForSVfit")
  process.DiTauPairsWHFit.config.event.resonances.Higgs.daughters.leg1.builder = process.nSVfitTauToHadBuilder
  process.DiTauPairsWHFit.config.event.resonances.Higgs.daughters.leg1.likelihoodFunctions = cms.VPSet(process.nSVfitTauLikelihoodPhaseSpace)
  process.DiTauPairsWHFit.config.event.resonances.Higgs.daughters.leg2.src = cms.InputTag("selectedPatTausForSVfit")
  process.DiTauPairsWHFit.config.event.resonances.Higgs.daughters.leg2.builder = process.nSVfitTauToHadBuilder
  process.DiTauPairsWHFit.config.event.resonances.Higgs.daughters.leg2.likelihoodFunctions = cms.VPSet(process.nSVfitTauLikelihoodPhaseSpace)
  process.DiTauPairsWHFit.config.event.resonances.Higgs.builder.polStates = cms.vstring('undefined')
  process.DiTauPairsWHFit.config.event.srcMEt = cms.InputTag('patMETsPF')
  #process.DiTauPairsWHFit.config.event.srcMEt = cms.InputTag('patType1CorrectedPFMet')

  #-------For WH tau+tau SVfit in Integration mode--------
  process.DiTauPairsWHInt = copy.deepcopy(process.nSVfitProducerByIntegrationWH)
  process.DiTauPairsWHInt.config.event.likelihoodFunctions[0].srcMEtCovMatrix = cms.InputTag('pfMEtSignCovMatrix')
  process.DiTauPairsWHInt.config.event.resonances.W.builder.daughters.chargedLepton.src = cms.InputTag("selectedPatMuonsForSVfit")
  process.DiTauPairsWHInt.config.event.resonances.W.daughters.chargedLepton.src = cms.InputTag("selectedPatMuonsForSVfit")
  process.DiTauPairsWHInt.config.event.resonances.Higgs.daughters.leg1.src = cms.InputTag("selectedPatTausForSVfit")
  process.DiTauPairsWHInt.config.event.resonances.Higgs.daughters.leg1.builder = process.nSVfitTauToHadBuilder
  process.DiTauPairsWHInt.config.event.resonances.Higgs.daughters.leg1.likelihoodFunctions = cms.VPSet(process.nSVfitTauLikelihoodPhaseSpace)
  process.DiTauPairsWHInt.config.event.resonances.Higgs.daughters.leg1.applySinThetaFactor = cms.bool(False)
  process.DiTauPairsWHInt.config.event.resonances.Higgs.daughters.leg2.src = cms.InputTag("selectedPatTausForSVfit")
  process.DiTauPairsWHInt.config.event.resonances.Higgs.daughters.leg2.builder = process.nSVfitTauToHadBuilder
  process.DiTauPairsWHInt.config.event.resonances.Higgs.daughters.leg2.likelihoodFunctions = cms.VPSet(process.nSVfitTauLikelihoodPhaseSpace)
  process.DiTauPairsWHInt.config.event.resonances.Higgs.daughters.leg2.applySinThetaFactor = cms.bool(False)
  process.DiTauPairsWHInt.config.event.resonances.Higgs.builder.polStates = cms.vstring('undefined')
  process.DiTauPairsWHInt.config.event.srcMEt = cms.InputTag('patMETsPF')
  #process.DiTauPairsWHInt.config.event.srcMEt = cms.InputTag('patType1CorrectedPFMet')

  process.SVND = cms.Sequence(
     process.selectedPatTausForSVfit +
     process.selectedPatMuonsForSVfit +
     process.ak5PFJetsNotOverlappingWithLeptons*process.pfCandsNotInJet*process.pfMEtSignCovMatrix*process.DiTauPairsWHInt
  )
  return process.SVND
