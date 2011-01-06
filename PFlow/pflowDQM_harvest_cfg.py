import FWCore.ParameterSet.Config as cms

process = cms.Process('PFlowHarvest')

fa = 'RelValQCD'
fb = 'FlatPt_15_3000_Fast'
fc = 'ParticleFlow'

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.dqmSaver.convention = 'Offline'
process.dqmSaver.workflow = '/%s/%s/%s' % (fa, fb, fc)
#process.dqmSaver.workflow = '/A/B/C'
process.dqmEnv.subSystemFolder = 'ParticleFlow'

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.222.2.6 $'),
    annotation = cms.untracked.string('step3_DT2_1 nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    fileMode = cms.untracked.string('FULLMERGE')
)
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:MEtoEDM_PFlow.root'),
    processingMode = cms.untracked.string('RunsAndLumis')
)
#process.load("Configuration.StandardSequences.EDMtoMEAtJobEnd_cff")

process.load("Validation.RecoParticleFlow.PFlowValidationClient_cff")
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')

process.p = cms.Path(process.EDMtoME*process.pfJetClient*process.pfMETClient*process.dqmEnv*process.dqmSaver)

