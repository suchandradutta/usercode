#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"

ClassImp(Event) 
ClassImp(GenEvent) 
ClassImp(Electron) 
ClassImp(GenParticle) 
ClassImp(GenJet) 
ClassImp(MET) 
ClassImp(GenMET) 
ClassImp(Tau) 
ClassImp(CaloJet)
ClassImp(Muon)
ClassImp(Jet)
ClassImp(SuperCluster)
ClassImp(Vertex)
ClassImp(Trigger)


Event::Event() :
  run(0), 
  event(0), 
  lumis(0), 
  bunch(0), 
  orbit(0), 
  time(-1), 
  isdata(false), 
  isPhysDeclared(false),
  isBPTX0(false),
  isBSCMinBias(false),
  isBSCBeamHalo(false),
  isPrimaryVertex(false),
  isBeamScraping(false),
  passHBHENoiseFilter(false) 
{
  nPU.clear();
  bunchCrossing.clear();
}

GenEvent::GenEvent() :
  processID(0),
  ptHat(-999)
{
  pdfWeights.clear();
}
Electron::Electron() :
  eta(-999),
  phi(-999),
  pt(-999),
  trackPt(-999),
  energy(-999),
  caloEnergy(-999),
  charge(-9),
  hoe(-999),
  sigmaEtaEta(-999),
  sigmaIEtaIEta(-999),
  deltaPhiTrkSC(-999),
  deltaEtaTrkSC(-999),
  classif(-1),
  e1x5overe5x5(-999),
  e2x5overe5x5(-99), 
  isoEcal03(-999),
  isoHcal03(-999),
  isoTrk03(-999),
  isoEcal04(-999),
  isoHcal04(-999),
  isoTrk04(-999),
  isoRel03(-999),
  isoRel04(-999),
  missingHits(-1),
  dist_vec(-999),
  dCotTheta(-999),
  scEn(-999),
  scEta(-999),
  scPhi(-999),
  scET(-999),
  scRawEnergy(-999),
  vtxDist3D(-999),
  vtxIndex(-1),
  vtxDistZ(-999),
  pfRelIso(-999), 
  scE1E9(-999),
  scS4S1(-999),
  sckOutOfTime(-999),
  scEcalIso(-999),
  scHEEPEcalIso(-999),
  scHEEPTrkIso(-999)
{}

GenParticle::GenParticle() :
  eta(-999),
  phi(-999),
  p(-999),
  px(-999),
  py(-999),
  pz(-999),
  pt(-999),
  energy(-999),
  pdgId(-999),
  vx(-999),
  vy(-999),
  vz(-999),
  numDaught(-1),
  status(-999),
  motherIndex(-999) {}

GenJet::GenJet() :
  eta(-999),
  phi(-999),
  p(-999),
  pt(-999),
  energy(-999),
  emf(-999),
  hadf(-999) {}

MET::MET() :
  met(-999),
  metphi(-999),
  sumet(-999),
  metuncorr(-999),
  metphiuncorr(-999),
  sumetuncorr(-999) {}

GenMET::GenMET() :
  met(-999),
  metphi(-999),
  sumet(-999) {}

Tau::Tau() :
  eta(-999),
  phi(-999),
  pt(-999),
  energy(-999),
  charge(-999),
  mass(-999),
  leadChargedParticlePt(-999),
  leadNeutralParticlePt(-999),
  leadParticlePt(-999),
  numChargedHadronsSignalCone(-1),
  numNeutralHadronsSignalCone(-1),
  numPhotonsSignalCone(-1),
  numParticlesSignalCone(-1),
  numChargedHadronsIsoCone(-1),
  numNeutralHadronsIsoCone(-1),
  numPhotonsIsoCone(-1),
  numParticlesIsoCone(-1),
  ptSumPFChargedHadronsIsoCone(-999),
  ptSumPhotonsIsoCone(-999),
  decayModeFinding(-1),
  looseIsolation(-1),
  mediumIsolation(-1),
  tightIsolation(-1),
  againstMuonLoose(-1),
  againstMuonTight(-1),
  againstElectronLoose(-1), 
  againstElectronMedium(-1), 
  againstElectronTight(-1), 
  pfElectronMVA(-999),
  jetPt(-999),
  jetEta(-999),
  jetPhi(-999),
  ecalStripSumEOverPLead(-999),
  bremsRecoveryEOverPLead(-999),
  hcal3x3OverPLead(-999),
  etaetaMoment(-999),
  phiphiMoment(-999),
  vx(-999), vy(-999), vz(-999),
  zvertex(-999), ltsipt(-999) 
  {}

CaloJet::CaloJet() :
  eta(-999),
  phi(-999),
  pt(-999),
  pt_raw(-999),
  energy(-999),
  energy_raw(-999),
  jecUnc(-999),
  resJEC(-999),
  overlaps(-1),
  partonFlavour(-1),
  emf(-999),
  resEmf(-999),
  hadf(-999),
  n90Hits(-1),
  fHPD(-999),
  fRBX(-999),
  sigmaEta(-999),
  sigmaPhi(-999),
  trackCountingHighEffBTag(-999),
  trackCountingHighPurBTag(-999),
  simpleSecondaryVertexHighEffBTag(-999),
  simpleSecondaryVertexHighPurBTag(-999),
  jetProbabilityBTag(-999),
  jetBProbabilityBTag(-999),
  passLooseID(-1),
  passTightID(-1) {}

Muon::Muon() :
  eta(-999),
  phi(-999),
  pt(-999),
  p(-999),
  energy(-999),
  charge(-999),
  trkD0(-999),
  trkD0Error(-999),
  trkDz(-999),
  trkDzError(-999),
  globalChi2(-999),
  trkIso(-999),
  ecalIso(-999),
  hcalIso(-999),
  hoIso(-999),
  relIso(-999),
  passID(-1),
  vtxDist3D(-999),
  vtxIndex(-1),
  vtxDistZ(-999),
  pixHits(-1),
  trkHits(-1),
  matches(-1),
  pfRelIso(-999),
  isTrackerMuon(-1) {}

Jet::Jet() :
  eta(-999),
  phi(-999),
  pt(-999),
  pt_raw(-999),
  energy(-999),
  energy_raw(-999),
  jecUnc(-999),
  resJEC(-999),
  partonFlavour(-1),
  chargedEmEnergyFraction(-999),
  chargedHadronEnergyFraction(-999),
  chargedMuEnergyFraction(-999),
  electronEnergyFraction(-999),
  muonEnergyFraction(-999),
  neutralEmEnergyFraction(-999),
  neutralHadronEnergyFraction(-999),
  photonEnergyFraction(-999),
  chargedHadronMultiplicity(-1),
  chargedMultiplicity(-1),
  electronMultiplicity(-1),
  muonMultiplicity(-1),
  neutralHadronMultiplicity(-1),
  neutralMultiplicity(-1),
  photonMultiplicity(-1),
  nConstituents(-1),
  trackCountingHighEffBTag(-999),
  trackCountingHighPurBTag(-999),
  simpleSecondaryVertexHighEffBTag(-999),
  simpleSecondaryVertexHighPurBTag(-999),
  jetProbabilityBTag(-999),
  jetBProbabilityBTag(-999),
  passLooseID(-1),
  passTightID(-1) {}

SuperCluster::SuperCluster() :
  eta(-999),
  phi(-999),
  pt(-999),
  rawEnergy(-999),
  clustersSize(-999),
  hoe(-999),
  sigmaIEtaIEta(-999),
  ecalIso(-999),
  e1e9(-999),
  s4s1(-999),
  kOutOfTime(-1),
  heepEcalIso(-999),
  heepTrkIso(-999),
  trackMatch(-1),
  dRTrack1(-999),
  track1Eta(-999),
  track1Phi(-999),
  track1Pt(-999),
  dRTrack2(-999),
  track2Eta(-999),
  track2Phi(-999),
  track2Pt(-999),
  scE1E9(-999),
  scS4S1(-999),
  sckOutOfTime(-1),
  scEcalIso(-999),
  scHEEPEcalIso(-999),
  scHEEPTrkIso(-999) {}

Vertex::Vertex() :
  x(-999),
  y(-999),
  z(-999),
  xErr(-999),
  yErr(-999),
  zErr(-999),
  rho(-999),
  chi2(-999),
  ndf(-1),
  ntracks(-1),
  ntracksw05(-1),
  isfake(true) {}

Trigger::Trigger() 
{
  l1physbits.clear();
  l1techbits.clear();
  hltbits.clear();
  hltresults.clear();
  hltprescales.clear();
}
