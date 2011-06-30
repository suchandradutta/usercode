#ifndef __PhysicsObjects_h
#define __PhysicsObjects_h

#include "vector"
#include "TObject.h"

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

  ClassDef(Event,1) 
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
  double trackPt;
  double energy;
  double caloEnergy;
     int charge;

  // ID variables
  double hoe;
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
  double pfRelIso;

  double scE1E9;
  double scS4S1;
  double sckOutOfTime;
  double scEcalIso;
  double scHEEPEcalIso;
  double scHEEPTrkIso;

  ClassDef(Electron,1) 
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
     int numDaught;
     int status;
     int motherIndex;

  ClassDef(GenParticle,1) 
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
     int decayModeFinding;
     int looseIsolation;
     int mediumIsolation;
     int tightIsolation;

  // discriminators against electrons/muons
     int againstMuonLoose;
     int againstMuonTight;
     int againstElectronLoose; 
     int againstElectronMedium; 
     int againstElectronTight; 
  double pfElectronMVA;

  // kinematic variables for PFJet associated to PFTau
  double jetPt;
  double jetEta;
  double jetPhi;

  double ecalStripSumEOverPLead;
  double bremsRecoveryEOverPLead;
  double hcal3x3OverPLead;

  double etaetaMoment;
  double phiphiMoment;

  double vx;
  double vy;
  double vz;

  double zvertex;
  double ltsipt;

  ClassDef(Tau,2) 
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

  ClassDef(CaloJet, 1)
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
     int matches;
  double pfRelIso;
     int isTrackerMuon;

  ClassDef(Muon, 1)
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

  ClassDef(Jet, 1)
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

  ClassDef(SuperCluster, 1)
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

  ClassDef(Vertex, 1)
};
class Trigger: public TObject {
public:
  Trigger();
  ~Trigger() {}

  std::vector<int> l1physbits;
  std::vector<int> l1techbits;
  std::vector<int> hltbits;
  std::vector<int> hltresults;
  std::vector<int> hltprescales;

  ClassDef(Trigger, 1)
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
#endif
