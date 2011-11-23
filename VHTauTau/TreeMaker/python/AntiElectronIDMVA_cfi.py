import FWCore.ParameterSet.Config as cms

antiElectronMVAIDcfg = cms.PSet(
  methodName = cms.string("BDT"),
  weightsX0BL = cms.FileInPath(
    "UserCode/Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_X_0BL_BDT.weights.xml"
  ),
  weights11BL = cms.FileInPath(
    "UserCode/Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_1_1BL_BDT.weights.xml"
  ),
  weights01BL = cms.FileInPath(
    "UserCode/Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_0_1BL_BDT.weights.xml"
  ),
  weightsX0EC = cms.FileInPath(
    "UserCode/Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_X_0EC_BDT.weights.xml"
  ),
  weights11EC = cms.FileInPath(
    "UserCode/Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_1_1EC_BDT.weights.xml"
  ),
  weights01EC = cms.FileInPath(
    "UserCode/Bianchi/Utilities/data/antiE_v4/TMVAClassification_v2_0_1EC_BDT.weights.xml"
  )
)
