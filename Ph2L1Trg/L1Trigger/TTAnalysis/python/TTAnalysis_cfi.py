import FWCore.ParameterSet.Config as cms

ttAnalysis = cms.EDAnalyzer("TTAnalysis",
    L1VertexInputTag = cms.InputTag("L1VertexFinder", "l1vertices"),
    L1TrackInputTag = cms.InputTag("l1tTTTracksFromTracklet", "Level1TTTracks"),
    L1EmulatedTrackInputTag = cms.InputTag("l1tTTTracksFromTrackletEmulation", "Level1TTTracks"),
    L1EGammaInputTag = cms.VInputTag(cms.InputTag("l1tEGammaClusterEmuProducer",""),
                                     cms.InputTag("l1tLayer1EG","L1EgEE")),
    L1TkEmuElectronInputTag = cms.VInputTag(cms.InputTag("l1tLayer1EG", "L1TkEleEB"),
                                            cms.InputTag("l1tLayer1EG", "L1TkEleEE")),
    L1TkEmuEmInputTag       = cms.VInputTag(cms.InputTag("l1tLayer1EG", "L1TkEmEB"),
                                            cms.InputTag("l1tLayer1EG", "L1TkEmEE")),
    L1TkVertexInputTag    = cms.InputTag("L1TkPrimaryVertex",""),
    L1EmuVertexInputTag   = cms.InputTag("l1tVertexFinderEmulator","l1verticesEmulation"),
    L1PFCandidateInputTag = cms.InputTag("l1tLayer1", "PF"),
    L1PuppiCandidateInputTag = cms.InputTag("l1tLayer1", "Puppi"),
    OfflineTrackInputTag  = cms.InputTag("generalTracks",""),
    OfflineVertexInputTag = cms.InputTag("goodOfflinePrimaryVertices"),
    BeamSpotInputTag      = cms.InputTag("offlineBeamSpot"),
    GenParticleInputTag   = cms.InputTag("genParticles",""),
    DebugFlag             = cms.int32(1),
    L1TrackFlag           = cms.bool(False),
    L1EmulatedTrackFlag   = cms.bool(False),
    L1EGammaFlag          = cms.bool(False),
    L1VertexFlag          = cms.bool(False),
    OfflineVertexFlag     = cms.bool(False),
    L1EmuVertexFlag       = cms.bool(False),
    L1TkVertexFlag        = cms.bool(False),
    L1TkEmuEmFlag         = cms.bool(False),
    L1TkEmuElectronFlag   = cms.bool(False),
    L1PFCandidateFlag     = cms.bool(False),
    L1PuppiCandidateFlag  = cms.bool(False),
    OfflineTrackFlag      = cms.bool(False),
    SimTrackFlag          = cms.bool(False),
    GenParticleFlag       = cms.bool(False)
)    
