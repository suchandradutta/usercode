import FWCore.ParameterSet.Config as cms

muonIdMVAcfg = cms.PSet(
  IdBarrel1Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml"),
  IdEndcap1Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml"),
  IdBarrel2Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml"),
  IdEndcap2Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml"),
  IdTrackerWeights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Tracker_V0_BDTG.weights.xml"),
  IdGlobalWeights  = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIDMVA_sixie-Global_V0_BDTG.weights.xml"),
  target           = cms.string("2011Data")
)
muonIsoMVAcfg = cms.PSet(
  IsoBarrel1Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml"),
  IsoEndcap1Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml"),
  IsoBarrel2Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml"),
  IsoEndcap2Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml"),
  IsoTrackerWeights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml"),
  IsoGlobalWeights  = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml")
)
muonIsoRingsRadMVAcfg = cms.PSet(
  IsoRingsRadBarrel1Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt5To10_V1_BDTG.weights.xml"),
  IsoRingsRadEndcap1Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt5To10_V1_BDTG.weights.xml"),
  IsoRingsRadBarrel2Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-BarrelPt10ToInf_V1_BDTG.weights.xml"),
  IsoRingsRadEndcap2Weights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-EndcapPt10ToInf_V1_BDTG.weights.xml"),
  IsoRingsRadTrackerWeights = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Tracker_V1_BDTG.weights.xml"),
  IsoRingsRadGlobalWeights  = cms.FileInPath("Muon/MuonAnalysisTools/data/MuonIsoMVA_sixie-Global_V1_BDTG.weights.xml")
)
