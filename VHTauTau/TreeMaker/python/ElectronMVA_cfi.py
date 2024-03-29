import FWCore.ParameterSet.Config as cms

electronIdMVAcfg = cms.PSet(
  IdCat1Weights = cms.FileInPath("EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat1.weights.xml"),
  IdCat2Weights = cms.FileInPath("EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat2.weights.xml"),
  IdCat3Weights = cms.FileInPath("EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat3.weights.xml"),
  IdCat4Weights = cms.FileInPath("EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat4.weights.xml"),
  IdCat5Weights = cms.FileInPath("EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat5.weights.xml"),
  IdCat6Weights = cms.FileInPath("EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat6.weights.xml")
)
electronIsoMVAcfg = cms.PSet(
  IsoBarrelPt1Weights = cms.FileInPath("UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_BarrelPt5To10.weights.xml"),
  IsoEndcapPt1Weights = cms.FileInPath("UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_EndcapPt5To10.weights.xml"),
  IsoBarrelPt2Weights = cms.FileInPath("UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_BarrelPt10ToInf.weights.xml"),
  IsoEndcapPt2Weights = cms.FileInPath("UserCode/sixie/EGamma/EGammaAnalysisTools/data/ElectronIso_BDTG_V0_EndcapPt10ToInf.weights.xml"),
  target              = cms.string("2011Data")
)
