#ifndef __VHTauTau_TreeMaker_PhysicsObjects_h
#define __VHTauTau_TreeMaker_PhysicsObjects_h

#include <vector>
#include <map>
#include <string>

#include "TObject.h"

namespace vhtm {
  class Event: public TObject {
  public:
    Event();
    ~Event() {}
  
    unsigned int run;
    unsigned int event;
    unsigned int lumis;
    unsigned int bunch;
    unsigned int orbit;
          double time;
            bool isdata;
  
    bool isPhysDeclared;
    bool isBPTX0;
    bool isBSCMinBias;
    bool isBSCBeamHalo;
    bool isPrimaryVertex;
    bool isBeamScraping;
    bool passHBHENoiseFilter;
  
    std::vector<int> nPU;
    std::vector<int> bunchCrossing;
    std::vector<int> trueNInt;
  
    ClassDef(Event,2) 
  };
  class GenEvent: public TObject {
  public:
    GenEvent();
    ~GenEvent() {}
  
    unsigned int processID;
    double ptHat;
    std::vector<double> pdfWeights;
  
    ClassDef(GenEvent,1) 
  };
  
  class Electron: public TObject {
  public:
    Electron();
    ~Electron() {}
    double eta;
    double phi;
    double pt;
      bool hasGsfTrack;
    double trackPt;
    double energy;
    double caloEnergy;
       int charge;
       int nValidHits;
     float simpleEleId60cIso;
     float simpleEleId70cIso;
     float simpleEleId80cIso;
     float simpleEleId85cIso;
     float simpleEleId90cIso;
     float simpleEleId95cIso;
  
    // ID variables
    double hoe;
    double eop;
    double sigmaEtaEta;
    double sigmaIEtaIEta;
    double deltaPhiTrkSC;
    double deltaEtaTrkSC;
       int classif;
    double e1x5overe5x5;
    double e2x5overe5x5; 
  
    // Iso variables
    double isoEcal03;
    double isoHcal03;
    double isoTrk03;
    double isoEcal04;
    double isoHcal04;
    double isoTrk04;
    double isoRel03;
    double isoRel04;
  
    // Conversion variables
       int missingHits;
    double dist_vec;
    double dCotTheta;
  
    // SC associated with electron
    double scEn;
    double scEta;
    double scPhi;
    double scET;
    double scRawEnergy;
  
    // Vertex association variables
    double vtxDist3D;
       int vtxIndex;
    double vtxDistZ;
    double relIso;
    double pfRelIso;
  
    double dB;
    double edB;
  
    double scE1E9;
    double scS4S1;
    double sckOutOfTime;
    double scEcalIso;
    double scHEEPEcalIso;
    double scHEEPTrkIso;
  
       int nBrems;
     float fbrem;
  
    double mva;
   
       int selbit;

    ClassDef(Electron, 4) 
  };
  class GenParticle: public TObject {
  public:
    GenParticle();
    ~GenParticle() {}
  
    double eta;
    double phi;
    double p;
    double px;
    double py;
    double pz;
    double pt;
    double energy;
       int pdgId;
    double vx;
    double vy;
    double vz;
       int status;
       int numDaught;
       int numMother;
       int motherIndex;
    std::vector<int> motherIndices;
    std::vector<int> daughtIndices;
  
    ClassDef(GenParticle,3) 
  };
  class GenJet: public TObject {
  public:
    GenJet();
    ~GenJet() {}
  
    double eta;
    double phi;
    double p;
    double pt;
    double energy;
    double emf;
    double hadf;
  
    ClassDef(GenJet,1) 
  };
  class MET: public TObject {
  public:
    MET();
    ~MET() {}
  
    double met;
    double metphi;
    double sumet;
    double metuncorr;
    double metphiuncorr;
    double sumetuncorr;
  
    ClassDef(MET,1) 
  };
  class Tau: public TObject {
  public:
    Tau();
    ~Tau() {}
  
    double eta;
    double phi;
    double pt;
    double energy;
       int charge;
    double mass;
  
    // Leading particle pT
    double leadChargedParticlePt;
    double leadNeutralParticlePt;
    double leadParticlePt;
  
    // Number of charged/neutral candidates and photons in different cones
       int numChargedHadronsSignalCone;
       int numNeutralHadronsSignalCone;
       int numPhotonsSignalCone;
       int numParticlesSignalCone;
       int numChargedHadronsIsoCone;
       int numNeutralHadronsIsoCone;
       int numPhotonsIsoCone;
       int numParticlesIsoCone;
    double ptSumPFChargedHadronsIsoCone;
    double ptSumPhotonsIsoCone;
    
    // tau id. discriminators
     float decayModeFinding;
     float looseIsolation;
     float mediumIsolation;
     float tightIsolation;
  
    // discriminators against electrons/muons
    float againstMuonLoose;
    float againstMuonTight;
    float againstElectronLoose; 
    float againstElectronMedium; 
    float againstElectronTight; 
    float pfElectronMVA;
    float againstElectronMVA;
  
    float byVLooseCombinedIsolationDeltaBetaCorr;
    float byLooseCombinedIsolationDeltaBetaCorr;
    float byMediumCombinedIsolationDeltaBetaCorr;
    float byTightCombinedIsolationDeltaBetaCorr;
    float byVLooseIsolationDeltaBetaCorr;
    float byLooseIsolationDeltaBetaCorr;
    float byMediumIsolationDeltaBetaCorr;
    float byTightIsolationDeltaBetaCorr;
  
    // kinematic variables for PFJet associated to PFTau
    double jetPt;
    double jetEta;
    double jetPhi;
  
    float maximumHCALPFClusterEt;
    float ecalStripSumEOverPLead;
    float bremsRecoveryEOverPLead;
    float hcalTotOverPLead;
    float hcalMaxOverPLead;
    float hcal3x3OverPLead;
  
    float etaetaMoment;
    float phiphiMoment;
    float etaphiMoment;
  
    double vx;
    double vy;
    double vz;
  
    double zvertex;
    double ltsipt;
  
    double mva;
    
       int selbit;
  
    ClassDef(Tau, 4) 
  };
  class CaloJet: public TObject {
  public:
    CaloJet();
    ~CaloJet() {}
    double eta;
    double phi;
    double pt;
    double pt_raw;
    double energy;
    double energy_raw;
    double jecUnc;
    double resJEC;
       int overlaps;
       int partonFlavour;
    double emf;
    double resEmf;
    double hadf;
       int n90Hits;
    double fHPD;
    double fRBX;
    double sigmaEta;
    double sigmaPhi;
    double trackCountingHighEffBTag;
    double trackCountingHighPurBTag;
    double simpleSecondaryVertexHighEffBTag;
    double simpleSecondaryVertexHighPurBTag;
    double jetProbabilityBTag;
    double jetBProbabilityBTag;
       int passLooseID;
       int passTightID;

       int selbit;
  
    ClassDef(CaloJet, 2)
  };
  class Muon: public TObject {
  public:
    Muon();
    ~Muon() {}
    double eta;
    double phi;
    double pt;
    double p;
    double energy;
       int charge;
    double trkD0;
    double trkD0Error;
    double trkDz;
    double trkDzError;
    double globalChi2;
    double trkIso;
    double ecalIso;
    double hcalIso;
    double hoIso;
    double relIso;
       int passID;
    double vtxDist3D;
       int vtxIndex;
    double vtxDistZ;
       int pixHits;
       int trkHits;
       int muoHits;
       int matches;
    double pfRelIso;
      bool isTrackerMuon;
    double dB;  // PV2D
    double edB;   
  
    // UW Recommendation
      bool isGlobalMuonPromptTight;
      bool isAllArbitrated;
       int nChambers;
       int nMatches;
       int nMatchedStations;
       unsigned int stationMask;
       unsigned int stationGapMaskDistance;
       unsigned int stationGapMaskPull;
     
       int selbit;

    ClassDef(Muon, 4)
  };
  class Jet: public TObject {
  public:
    Jet();
    ~Jet() {}
    double eta;
    double phi;
    double pt;
    double pt_raw;
    double energy;
    double energy_raw;
    double jecUnc;
    double resJEC;
       int partonFlavour;
    double chargedEmEnergyFraction;
    double chargedHadronEnergyFraction;
    double chargedMuEnergyFraction;
    double electronEnergyFraction;
    double muonEnergyFraction;
    double neutralEmEnergyFraction;
    double neutralHadronEnergyFraction;
    double photonEnergyFraction;
       int chargedHadronMultiplicity;
       int chargedMultiplicity;
       int electronMultiplicity;
       int muonMultiplicity;
       int neutralHadronMultiplicity;
       int neutralMultiplicity;
       int photonMultiplicity;
       int nConstituents;
    double trackCountingHighEffBTag;
    double trackCountingHighPurBTag;
    double simpleSecondaryVertexHighEffBTag;
    double simpleSecondaryVertexHighPurBTag;
    double jetProbabilityBTag;
    double jetBProbabilityBTag;
       int passLooseID;
       int passTightID;
  
       int selbit;

    ClassDef(Jet, 2)
  };
  class SuperCluster: public TObject {
  public: 
    SuperCluster();
    ~SuperCluster() {}
  
    double eta;
    double phi;
    double pt;
    double rawEnergy;
       int clustersSize;
    double hoe;
    double sigmaIEtaIEta;
    double ecalIso;
    double e1e9;
    double s4s1;
       int kOutOfTime;
    double heepEcalIso;
    double heepTrkIso;
       int trackMatch;
    double dRTrack1;
    double track1Eta;
    double track1Phi;
    double track1Pt;
    double dRTrack2;
    double track2Eta;
    double track2Phi;
    double track2Pt;
    double scE1E9;
    double scS4S1;
       int sckOutOfTime;
    double scEcalIso;
    double scHEEPEcalIso;
    double scHEEPTrkIso;
  
       int selbit;

    ClassDef(SuperCluster, 2)
  };
  class Vertex: public TObject {
  public:
    Vertex();
    ~Vertex() {}
  
    double x;
    double y;
    double z;
    double xErr;
    double yErr;
    double zErr;
    double rho;
    double chi2;
    double ndf;
       int ntracks;
       int ntracksw05;
      bool isfake;
      bool isvalid;
    double sumPt;
  
       int selbit;

    ClassDef(Vertex, 3)
  };
  class Trigger: public TObject {
  public:
    Trigger();
    ~Trigger() {}
  
    std::vector<int> l1physbits;
    std::vector<int> l1techbits;
    std::vector<std::string> hltpaths;
    std::vector<int> hltresults;
    std::vector<int> hltprescales;
  
    ClassDef(Trigger, 2)
  };
  class GenMET: public TObject {
  public:
    GenMET();
    ~GenMET() {}
  
    double met;
    double metphi;
    double sumet;
  
    ClassDef(GenMET, 1)
  };
  class Track: public TObject {
  public:
    Track();
    ~Track() {}
  
    double eta;
    double etaError;
    double theta;
    double thetaError;
    double phi;
    double phiError;
    double p;
    double pt;
    double ptError;
    double qoverp;
    double qoverpError;
  
       int nValidHits;
       int nLostHits;
    double validFraction;
       int nValidTrackerHits;
       int nValidPixelHits;
       int nValidStripHits;
       int trackerLayersWithMeasurement;
       int pixelLayersWithMeasurement;
       int stripLayersWithMeasurement;
  
       double dxy;
       double dxyError;
       double dz;
       double dzError;
  
       double chi2;
          int ndof;
       double vx;
       double vy;
       double vz;
  
          int selbit; 
         
    ClassDef(Track, 2)
  };
  class Photon : public TObject {
  public:
    Photon();
    ~Photon() {}
  
    double et;
    double eta;
    double phi;
    double energy;
    double theta; 
    double vx;
    double vy;
    double vz;
  
    double scEnergy;
    double scEta;
    double scPhi;
    double scSize;
    double scEtaWidth;
    double scPhiWidth;
    double scEt;
    double scRawEnergy;
    double scx;
    double scy;
    double scz; 
    double isoEcalRecHit03;
    double isoHcalRecHit03;
    double isoSolidTrkCone03;
    double isoHollowTrkCone03;
       int nTrkSolidCone03;
       int nTrkHollowCone03;
  
    double isoEcalRecHit04;
    double isoHcalRecHit04;
    double isoSolidTrkCone04;
    double isoHollowTrkCone04;
       int nTrkSolidCone04;
       int nTrkHollowCone04;
  
      bool isEB;
      bool isEE; 
      bool isEBGap;
      bool isEEGap;
      bool isEBEEGap;
  
      bool hasPixelSeed;
    double ecalIso;
    double hcalIso;
    double trackIso;
    double chargedHadIso;
    double neutralHadIso;
    double photonIso;
  
    double r9;
    double hoe;
    double sigmaEtaEta;
    double sigmaIEtaIEta;
    double e1x5;
    double e2x5; 
    double e3x3;
    double e5x5;
    double r1x5;
    double r2x5;
    double maxEnergyXtal;
  
      bool hasConversionTracks;
       int nTracks;
      bool isConverted;
    double pairInvMass;
    double pairCotThetaSeparation;
    double pairPx;
    double pairPy;
    double pairPz;
    double conv_vx;
    double conv_vy;
    double conv_vz;
    double eovp;
    double zpv;
    double distOfMinApproach;
    double dPhiTracksAtVtx;
    double dPhiTracksAtEcal;
    double dEtaTracksAtEcal;  

       int selbit;
  
    ClassDef(Photon, 2)
  };
  class TriggerObject: public TObject {
  public:
    TriggerObject();
    ~TriggerObject() {}
  
    double energy;
    double pt;
    double eta;
    double phi;
  
    std::map< std::string, unsigned int > pathList;
  
    ClassDef(TriggerObject, 1)
  };
  class CommonVertex: public TObject {
  public:
    CommonVertex();
    ~CommonVertex() {}
  
    double chi2;
    int ndof;
    std::string label;
  
    ClassDef(CommonVertex, 1)
  };
}
#endif
