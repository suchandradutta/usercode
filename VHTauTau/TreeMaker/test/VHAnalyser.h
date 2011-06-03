//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat May 28 17:31:50 2011 by ROOT version 5.27/06b
// from TTree vhtree/
// found on file: MC_RelValTTbar.root
//////////////////////////////////////////////////////////

#ifndef VHAnalyser_h
#define VHAnalyser_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

   const Int_t kMaxEvent = 1;
   const Int_t kMaxVertex = 2;
   const Int_t kMaxCaloJet = 17;
   const Int_t kMaxJet = 1;
   const Int_t kMaxElectron = 5;
   const Int_t kMaxMET = 1;
   const Int_t kMaxMuon = 2;
   const Int_t kMaxTau = 1;
   const Int_t kMaxGenParticle = 500;
   const Int_t kMaxGenJet = 34;
   const Int_t kMaxGenMET = 1;

class VHAnalyser {
public :

   // User Code 
   void bookHistograms();
   void saveHistograms();
   
   Bool_t          bookedHistos;
   std::string     fileName;
   TFile*          outputFile;

   TH1*            hNVertex;


   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Event_;
   UInt_t          Event_fUniqueID[kMaxEvent];   //[Event_]
   UInt_t          Event_fBits[kMaxEvent];   //[Event_]
   UInt_t          Event_run[kMaxEvent];   //[Event_]
   UInt_t          Event_event[kMaxEvent];   //[Event_]
   UInt_t          Event_lumis[kMaxEvent];   //[Event_]
   UInt_t          Event_bunch[kMaxEvent];   //[Event_]
   UInt_t          Event_orbit[kMaxEvent];   //[Event_]
   Double_t        Event_time[kMaxEvent];   //[Event_]
   Bool_t          Event_isdata[kMaxEvent];   //[Event_]
   Int_t           Vertex_;
   UInt_t          Vertex_fUniqueID[kMaxVertex];   //[Vertex_]
   UInt_t          Vertex_fBits[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_x[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_y[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_z[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_xErr[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_yErr[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_zErr[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_rho[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_chi2[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_ndf[kMaxVertex];   //[Vertex_]
   Int_t           Vertex_ntracks[kMaxVertex];   //[Vertex_]
   Int_t           Vertex_ntracksw05[kMaxVertex];   //[Vertex_]
   Bool_t          Vertex_isfake[kMaxVertex];   //[Vertex_]
   Int_t           nVertex;
   Int_t           CaloJet_;
   UInt_t          CaloJet_fUniqueID[kMaxCaloJet];   //[CaloJet_]
   UInt_t          CaloJet_fBits[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_eta[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_phi[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_pt[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_pt_raw[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_energy[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_energy_raw[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_jecUnc[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_resJEC[kMaxCaloJet];   //[CaloJet_]
   Int_t           CaloJet_overlaps[kMaxCaloJet];   //[CaloJet_]
   Int_t           CaloJet_partonFlavour[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_emf[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_resEmf[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_hadf[kMaxCaloJet];   //[CaloJet_]
   Int_t           CaloJet_n90Hits[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_fHPD[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_fRBX[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_sigmaEta[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_sigmaPhi[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_trackCountingHighEffBTag[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_trackCountingHighPurBTag[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_simpleSecondaryVertexHighEffBTag[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_simpleSecondaryVertexHighPurBTag[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_jetProbabilityBTag[kMaxCaloJet];   //[CaloJet_]
   Double_t        CaloJet_jetBProbabilityBTag[kMaxCaloJet];   //[CaloJet_]
   Int_t           CaloJet_passLooseID[kMaxCaloJet];   //[CaloJet_]
   Int_t           CaloJet_passTightID[kMaxCaloJet];   //[CaloJet_]
   Int_t           nCaloJet;
   Int_t           Jet_;
   UInt_t          Jet_fUniqueID[kMaxJet];   //[Jet_]
   UInt_t          Jet_fBits[kMaxJet];   //[Jet_]
   Double_t        Jet_eta[kMaxJet];   //[Jet_]
   Double_t        Jet_phi[kMaxJet];   //[Jet_]
   Double_t        Jet_pt[kMaxJet];   //[Jet_]
   Double_t        Jet_pt_raw[kMaxJet];   //[Jet_]
   Double_t        Jet_energy[kMaxJet];   //[Jet_]
   Double_t        Jet_energy_raw[kMaxJet];   //[Jet_]
   Double_t        Jet_jecUnc[kMaxJet];   //[Jet_]
   Double_t        Jet_resJEC[kMaxJet];   //[Jet_]
   Int_t           Jet_partonFlavour[kMaxJet];   //[Jet_]
   Double_t        Jet_chargedEmEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_chargedHadronEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_chargedMuEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_electronEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_muonEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_neutralEmEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_neutralHadronEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_photonEnergyFraction[kMaxJet];   //[Jet_]
   Int_t           Jet_chargedHadronMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_chargedMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_electronMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_muonMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_neutralHadronMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_neutralMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_photonMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_nConstituents[kMaxJet];   //[Jet_]
   Double_t        Jet_trackCountingHighEffBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_trackCountingHighPurBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_simpleSecondaryVertexHighEffBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_simpleSecondaryVertexHighPurBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_jetProbabilityBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_jetBProbabilityBTag[kMaxJet];   //[Jet_]
   Int_t           Jet_passLooseID[kMaxJet];   //[Jet_]
   Int_t           Jet_passTightID[kMaxJet];   //[Jet_]
   Int_t           nJet;
   Int_t           Electron_;
   UInt_t          Electron_fUniqueID[kMaxElectron];   //[Electron_]
   UInt_t          Electron_fBits[kMaxElectron];   //[Electron_]
   Double_t        Electron_eta[kMaxElectron];   //[Electron_]
   Double_t        Electron_phi[kMaxElectron];   //[Electron_]
   Double_t        Electron_pt[kMaxElectron];   //[Electron_]
   Double_t        Electron_trackPt[kMaxElectron];   //[Electron_]
   Double_t        Electron_energy[kMaxElectron];   //[Electron_]
   Double_t        Electron_caloEnergy[kMaxElectron];   //[Electron_]
   Int_t           Electron_charge[kMaxElectron];   //[Electron_]
   Double_t        Electron_hoe[kMaxElectron];   //[Electron_]
   Double_t        Electron_sigmaEtaEta[kMaxElectron];   //[Electron_]
   Double_t        Electron_sigmaIEtaIEta[kMaxElectron];   //[Electron_]
   Double_t        Electron_deltaPhiTrkSC[kMaxElectron];   //[Electron_]
   Double_t        Electron_deltaEtaTrkSC[kMaxElectron];   //[Electron_]
   Int_t           Electron_classif[kMaxElectron];   //[Electron_]
   Double_t        Electron_e1x5overe5x5[kMaxElectron];   //[Electron_]
   Double_t        Electron_e2x5overe5x5[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoEcal03[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoHcal03[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoTrk03[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoEcal04[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoHcal04[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoTrk04[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoRel03[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoRel04[kMaxElectron];   //[Electron_]
   Int_t           Electron_missingHits[kMaxElectron];   //[Electron_]
   Double_t        Electron_dist_vec[kMaxElectron];   //[Electron_]
   Double_t        Electron_dCotTheta[kMaxElectron];   //[Electron_]
   Double_t        Electron_scEn[kMaxElectron];   //[Electron_]
   Double_t        Electron_scEta[kMaxElectron];   //[Electron_]
   Double_t        Electron_scPhi[kMaxElectron];   //[Electron_]
   Double_t        Electron_scET[kMaxElectron];   //[Electron_]
   Double_t        Electron_scRawEnergy[kMaxElectron];   //[Electron_]
   Double_t        Electron_vtxDist3D[kMaxElectron];   //[Electron_]
   Int_t           Electron_vtxIndex[kMaxElectron];   //[Electron_]
   Double_t        Electron_vtxDistZ[kMaxElectron];   //[Electron_]
   Double_t        Electron_pfRelIso[kMaxElectron];   //[Electron_]
   Int_t           nElectron;
   Int_t           MET_;
   UInt_t          MET_fUniqueID[kMaxMET];   //[MET_]
   UInt_t          MET_fBits[kMaxMET];   //[MET_]
   Double_t        MET_met[kMaxMET];   //[MET_]
   Double_t        MET_metphi[kMaxMET];   //[MET_]
   Double_t        MET_sumet[kMaxMET];   //[MET_]
   Double_t        MET_metuncorr[kMaxMET];   //[MET_]
   Double_t        MET_metphiuncorr[kMaxMET];   //[MET_]
   Double_t        MET_sumetuncorr[kMaxMET];   //[MET_]
   Int_t           nMET;
   Int_t           Muon_;
   UInt_t          Muon_fUniqueID[kMaxMuon];   //[Muon_]
   UInt_t          Muon_fBits[kMaxMuon];   //[Muon_]
   Double_t        Muon_eta[kMaxMuon];   //[Muon_]
   Double_t        Muon_phi[kMaxMuon];   //[Muon_]
   Double_t        Muon_pt[kMaxMuon];   //[Muon_]
   Double_t        Muon_p[kMaxMuon];   //[Muon_]
   Double_t        Muon_energy[kMaxMuon];   //[Muon_]
   Int_t           Muon_charge[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkD0[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkD0Error[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkDz[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkDzError[kMaxMuon];   //[Muon_]
   Double_t        Muon_globalChi2[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkIso[kMaxMuon];   //[Muon_]
   Double_t        Muon_ecalIso[kMaxMuon];   //[Muon_]
   Double_t        Muon_hcalIso[kMaxMuon];   //[Muon_]
   Double_t        Muon_hoIso[kMaxMuon];   //[Muon_]
   Double_t        Muon_relIso[kMaxMuon];   //[Muon_]
   Int_t           Muon_passID[kMaxMuon];   //[Muon_]
   Double_t        Muon_vtxDist3D[kMaxMuon];   //[Muon_]
   Int_t           Muon_vtxIndex[kMaxMuon];   //[Muon_]
   Double_t        Muon_vtxDistZ[kMaxMuon];   //[Muon_]
   Int_t           Muon_pixHits[kMaxMuon];   //[Muon_]
   Int_t           Muon_trkHits[kMaxMuon];   //[Muon_]
   Int_t           Muon_matches[kMaxMuon];   //[Muon_]
   Double_t        Muon_pfRelIso[kMaxMuon];   //[Muon_]
   Int_t           Muon_isTrackerMuon[kMaxMuon];   //[Muon_]
   Int_t           nMuon;
   Int_t           Tau_;
   UInt_t          Tau_fUniqueID[kMaxTau];   //[Tau_]
   UInt_t          Tau_fBits[kMaxTau];   //[Tau_]
   Double_t        Tau_eta[kMaxTau];   //[Tau_]
   Double_t        Tau_phi[kMaxTau];   //[Tau_]
   Double_t        Tau_pt[kMaxTau];   //[Tau_]
   Double_t        Tau_energy[kMaxTau];   //[Tau_]
   Int_t           Tau_charge[kMaxTau];   //[Tau_]
   Double_t        Tau_mass[kMaxTau];   //[Tau_]
   Double_t        Tau_leadChargedParticlePt[kMaxTau];   //[Tau_]
   Double_t        Tau_leadNeutralParticlePt[kMaxTau];   //[Tau_]
   Double_t        Tau_leadParticlePt[kMaxTau];   //[Tau_]
   Int_t           Tau_numChargedHadronsSignalCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numNeutralHadronsSignalCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numPhotonsSignalCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numParticlesSignalCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numChargedHadronsIsoCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numNeutralHadronsIsoCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numPhotonsIsoCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numParticlesIsoCone[kMaxTau];   //[Tau_]
   Double_t        Tau_ptSumPFChargedHadronsIsoCone[kMaxTau];   //[Tau_]
   Double_t        Tau_ptSumPhotonsIsoCone[kMaxTau];   //[Tau_]
   Int_t           Tau_decayModeFinding[kMaxTau];   //[Tau_]
   Int_t           Tau_looseIsolation[kMaxTau];   //[Tau_]
   Int_t           Tau_mediumIsolation[kMaxTau];   //[Tau_]
   Int_t           Tau_tightIsolation[kMaxTau];   //[Tau_]
   Int_t           Tau_againstMuonLoose[kMaxTau];   //[Tau_]
   Int_t           Tau_againstMuonTight[kMaxTau];   //[Tau_]
   Int_t           Tau_againstElectronLoose[kMaxTau];   //[Tau_]
   Int_t           Tau_againstElectronMedium[kMaxTau];   //[Tau_]
   Int_t           Tau_againstElectronTight[kMaxTau];   //[Tau_]
   Double_t        Tau_pfElectronMVA[kMaxTau];   //[Tau_]
   Double_t        Tau_jetPt[kMaxTau];   //[Tau_]
   Double_t        Tau_jetEta[kMaxTau];   //[Tau_]
   Double_t        Tau_jetPhi[kMaxTau];   //[Tau_]
   Int_t           nTau;
   Int_t           GenParticle_;
   UInt_t          GenParticle_fUniqueID[kMaxGenParticle];   //[GenParticle_]
   UInt_t          GenParticle_fBits[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_eta[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_phi[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_p[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_px[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_py[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_pz[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_pt[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_energy[kMaxGenParticle];   //[GenParticle_]
   Int_t           GenParticle_pdgId[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_vx[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_vy[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_vz[kMaxGenParticle];   //[GenParticle_]
   Int_t           GenParticle_numDaught[kMaxGenParticle];   //[GenParticle_]
   Int_t           GenParticle_status[kMaxGenParticle];   //[GenParticle_]
   Int_t           GenParticle_motherIndex[kMaxGenParticle];   //[GenParticle_]
   Int_t           nGenParticle;
   Int_t           GenJet_;
   UInt_t          GenJet_fUniqueID[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_fBits[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_eta[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_phi[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_p[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_pt[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_energy[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_emf[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_hadf[kMaxGenJet];   //[GenJet_]
   Int_t           nGenJet;
   Int_t           GenMET_;
   UInt_t          GenMET_fUniqueID[kMaxGenMET];   //[GenMET_]
   UInt_t          GenMET_fBits[kMaxGenMET];   //[GenMET_]
   Double_t        GenMET_met[kMaxGenMET];   //[GenMET_]
   Double_t        GenMET_metphi[kMaxGenMET];   //[GenMET_]
   Double_t        GenMET_sumet[kMaxGenMET];   //[GenMET_]
   Int_t           nGenMET;

   // List of branches
   TBranch        *b_Event_;   //!
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_run;   //!
   TBranch        *b_Event_event;   //!
   TBranch        *b_Event_lumis;   //!
   TBranch        *b_Event_bunch;   //!
   TBranch        *b_Event_orbit;   //!
   TBranch        *b_Event_time;   //!
   TBranch        *b_Event_isdata;   //!
   TBranch        *b_Vertex_;   //!
   TBranch        *b_Vertex_fUniqueID;   //!
   TBranch        *b_Vertex_fBits;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_Vertex_xErr;   //!
   TBranch        *b_Vertex_yErr;   //!
   TBranch        *b_Vertex_zErr;   //!
   TBranch        *b_Vertex_rho;   //!
   TBranch        *b_Vertex_chi2;   //!
   TBranch        *b_Vertex_ndf;   //!
   TBranch        *b_Vertex_ntracks;   //!
   TBranch        *b_Vertex_ntracksw05;   //!
   TBranch        *b_Vertex_isfake;   //!
   TBranch        *b_fnVertex;   //!
   TBranch        *b_CaloJet_;   //!
   TBranch        *b_CaloJet_fUniqueID;   //!
   TBranch        *b_CaloJet_fBits;   //!
   TBranch        *b_CaloJet_eta;   //!
   TBranch        *b_CaloJet_phi;   //!
   TBranch        *b_CaloJet_pt;   //!
   TBranch        *b_CaloJet_pt_raw;   //!
   TBranch        *b_CaloJet_energy;   //!
   TBranch        *b_CaloJet_energy_raw;   //!
   TBranch        *b_CaloJet_jecUnc;   //!
   TBranch        *b_CaloJet_resJEC;   //!
   TBranch        *b_CaloJet_overlaps;   //!
   TBranch        *b_CaloJet_partonFlavour;   //!
   TBranch        *b_CaloJet_emf;   //!
   TBranch        *b_CaloJet_resEmf;   //!
   TBranch        *b_CaloJet_hadf;   //!
   TBranch        *b_CaloJet_n90Hits;   //!
   TBranch        *b_CaloJet_fHPD;   //!
   TBranch        *b_CaloJet_fRBX;   //!
   TBranch        *b_CaloJet_sigmaEta;   //!
   TBranch        *b_CaloJet_sigmaPhi;   //!
   TBranch        *b_CaloJet_trackCountingHighEffBTag;   //!
   TBranch        *b_CaloJet_trackCountingHighPurBTag;   //!
   TBranch        *b_CaloJet_simpleSecondaryVertexHighEffBTag;   //!
   TBranch        *b_CaloJet_simpleSecondaryVertexHighPurBTag;   //!
   TBranch        *b_CaloJet_jetProbabilityBTag;   //!
   TBranch        *b_CaloJet_jetBProbabilityBTag;   //!
   TBranch        *b_CaloJet_passLooseID;   //!
   TBranch        *b_CaloJet_passTightID;   //!
   TBranch        *b_fnCaloJet;   //!
   TBranch        *b_Jet_;   //!
   TBranch        *b_Jet_fUniqueID;   //!
   TBranch        *b_Jet_fBits;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_pt_raw;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_energy_raw;   //!
   TBranch        *b_Jet_jecUnc;   //!
   TBranch        *b_Jet_resJEC;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_Jet_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_chargedMuEnergyFraction;   //!
   TBranch        *b_Jet_electronEnergyFraction;   //!
   TBranch        *b_Jet_muonEnergyFraction;   //!
   TBranch        *b_Jet_neutralEmEnergyFraction;   //!
   TBranch        *b_Jet_neutralHadronEnergyFraction;   //!
   TBranch        *b_Jet_photonEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronMultiplicity;   //!
   TBranch        *b_Jet_chargedMultiplicity;   //!
   TBranch        *b_Jet_electronMultiplicity;   //!
   TBranch        *b_Jet_muonMultiplicity;   //!
   TBranch        *b_Jet_neutralHadronMultiplicity;   //!
   TBranch        *b_Jet_neutralMultiplicity;   //!
   TBranch        *b_Jet_photonMultiplicity;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_trackCountingHighEffBTag;   //!
   TBranch        *b_Jet_trackCountingHighPurBTag;   //!
   TBranch        *b_Jet_simpleSecondaryVertexHighEffBTag;   //!
   TBranch        *b_Jet_simpleSecondaryVertexHighPurBTag;   //!
   TBranch        *b_Jet_jetProbabilityBTag;   //!
   TBranch        *b_Jet_jetBProbabilityBTag;   //!
   TBranch        *b_Jet_passLooseID;   //!
   TBranch        *b_Jet_passTightID;   //!
   TBranch        *b_fnJet;   //!
   TBranch        *b_Electron_;   //!
   TBranch        *b_Electron_fUniqueID;   //!
   TBranch        *b_Electron_fBits;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_trackPt;   //!
   TBranch        *b_Electron_energy;   //!
   TBranch        *b_Electron_caloEnergy;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_sigmaEtaEta;   //!
   TBranch        *b_Electron_sigmaIEtaIEta;   //!
   TBranch        *b_Electron_deltaPhiTrkSC;   //!
   TBranch        *b_Electron_deltaEtaTrkSC;   //!
   TBranch        *b_Electron_classif;   //!
   TBranch        *b_Electron_e1x5overe5x5;   //!
   TBranch        *b_Electron_e2x5overe5x5;   //!
   TBranch        *b_Electron_isoEcal03;   //!
   TBranch        *b_Electron_isoHcal03;   //!
   TBranch        *b_Electron_isoTrk03;   //!
   TBranch        *b_Electron_isoEcal04;   //!
   TBranch        *b_Electron_isoHcal04;   //!
   TBranch        *b_Electron_isoTrk04;   //!
   TBranch        *b_Electron_isoRel03;   //!
   TBranch        *b_Electron_isoRel04;   //!
   TBranch        *b_Electron_missingHits;   //!
   TBranch        *b_Electron_dist_vec;   //!
   TBranch        *b_Electron_dCotTheta;   //!
   TBranch        *b_Electron_scEn;   //!
   TBranch        *b_Electron_scEta;   //!
   TBranch        *b_Electron_scPhi;   //!
   TBranch        *b_Electron_scET;   //!
   TBranch        *b_Electron_scRawEnergy;   //!
   TBranch        *b_Electron_vtxDist3D;   //!
   TBranch        *b_Electron_vtxIndex;   //!
   TBranch        *b_Electron_vtxDistZ;   //!
   TBranch        *b_Electron_pfRelIso;   //!
   TBranch        *b_fnElectron;   //!
   TBranch        *b_MET_;   //!
   TBranch        *b_MET_fUniqueID;   //!
   TBranch        *b_MET_fBits;   //!
   TBranch        *b_MET_met;   //!
   TBranch        *b_MET_metphi;   //!
   TBranch        *b_MET_sumet;   //!
   TBranch        *b_MET_metuncorr;   //!
   TBranch        *b_MET_metphiuncorr;   //!
   TBranch        *b_MET_sumetuncorr;   //!
   TBranch        *b_fnMET;   //!
   TBranch        *b_Muon_;   //!
   TBranch        *b_Muon_fUniqueID;   //!
   TBranch        *b_Muon_fBits;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_p;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_trkD0;   //!
   TBranch        *b_Muon_trkD0Error;   //!
   TBranch        *b_Muon_trkDz;   //!
   TBranch        *b_Muon_trkDzError;   //!
   TBranch        *b_Muon_globalChi2;   //!
   TBranch        *b_Muon_trkIso;   //!
   TBranch        *b_Muon_ecalIso;   //!
   TBranch        *b_Muon_hcalIso;   //!
   TBranch        *b_Muon_hoIso;   //!
   TBranch        *b_Muon_relIso;   //!
   TBranch        *b_Muon_passID;   //!
   TBranch        *b_Muon_vtxDist3D;   //!
   TBranch        *b_Muon_vtxIndex;   //!
   TBranch        *b_Muon_vtxDistZ;   //!
   TBranch        *b_Muon_pixHits;   //!
   TBranch        *b_Muon_trkHits;   //!
   TBranch        *b_Muon_matches;   //!
   TBranch        *b_Muon_pfRelIso;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_fnMuon;   //!
   TBranch        *b_Tau_;   //!
   TBranch        *b_Tau_fUniqueID;   //!
   TBranch        *b_Tau_fBits;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_energy;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_mass;   //!
   TBranch        *b_Tau_leadChargedParticlePt;   //!
   TBranch        *b_Tau_leadNeutralParticlePt;   //!
   TBranch        *b_Tau_leadParticlePt;   //!
   TBranch        *b_Tau_numChargedHadronsSignalCone;   //!
   TBranch        *b_Tau_numNeutralHadronsSignalCone;   //!
   TBranch        *b_Tau_numPhotonsSignalCone;   //!
   TBranch        *b_Tau_numParticlesSignalCone;   //!
   TBranch        *b_Tau_numChargedHadronsIsoCone;   //!
   TBranch        *b_Tau_numNeutralHadronsIsoCone;   //!
   TBranch        *b_Tau_numPhotonsIsoCone;   //!
   TBranch        *b_Tau_numParticlesIsoCone;   //!
   TBranch        *b_Tau_ptSumPFChargedHadronsIsoCone;   //!
   TBranch        *b_Tau_ptSumPhotonsIsoCone;   //!
   TBranch        *b_Tau_decayModeFinding;   //!
   TBranch        *b_Tau_looseIsolation;   //!
   TBranch        *b_Tau_mediumIsolation;   //!
   TBranch        *b_Tau_tightIsolation;   //!
   TBranch        *b_Tau_againstMuonLoose;   //!
   TBranch        *b_Tau_againstMuonTight;   //!
   TBranch        *b_Tau_againstElectronLoose;   //!
   TBranch        *b_Tau_againstElectronMedium;   //!
   TBranch        *b_Tau_againstElectronTight;   //!
   TBranch        *b_Tau_pfElectronMVA;   //!
   TBranch        *b_Tau_jetPt;   //!
   TBranch        *b_Tau_jetEta;   //!
   TBranch        *b_Tau_jetPhi;   //!
   TBranch        *b_fnTau;   //!
   TBranch        *b_GenParticle_;   //!
   TBranch        *b_GenParticle_fUniqueID;   //!
   TBranch        *b_GenParticle_fBits;   //!
   TBranch        *b_GenParticle_eta;   //!
   TBranch        *b_GenParticle_phi;   //!
   TBranch        *b_GenParticle_p;   //!
   TBranch        *b_GenParticle_px;   //!
   TBranch        *b_GenParticle_py;   //!
   TBranch        *b_GenParticle_pz;   //!
   TBranch        *b_GenParticle_pt;   //!
   TBranch        *b_GenParticle_energy;   //!
   TBranch        *b_GenParticle_pdgId;   //!
   TBranch        *b_GenParticle_vx;   //!
   TBranch        *b_GenParticle_vy;   //!
   TBranch        *b_GenParticle_vz;   //!
   TBranch        *b_GenParticle_numDaught;   //!
   TBranch        *b_GenParticle_status;   //!
   TBranch        *b_GenParticle_motherIndex;   //!
   TBranch        *b_fnGenParticle;   //!
   TBranch        *b_GenJet_;   //!
   TBranch        *b_GenJet_fUniqueID;   //!
   TBranch        *b_GenJet_fBits;   //!
   TBranch        *b_GenJet_eta;   //!
   TBranch        *b_GenJet_phi;   //!
   TBranch        *b_GenJet_p;   //!
   TBranch        *b_GenJet_pt;   //!
   TBranch        *b_GenJet_energy;   //!
   TBranch        *b_GenJet_emf;   //!
   TBranch        *b_GenJet_hadf;   //!
   TBranch        *b_fnGenJet;   //!
   TBranch        *b_GenMET_;   //!
   TBranch        *b_GenMET_fUniqueID;   //!
   TBranch        *b_GenMET_fBits;   //!
   TBranch        *b_GenMET_met;   //!
   TBranch        *b_GenMET_metphi;   //!
   TBranch        *b_GenMET_sumet;   //!
   TBranch        *b_fnGenMET;   //!

  VHAnalyser(TTree *tree, std::string fname);
   virtual ~VHAnalyser();

   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef VHAnalyser_cxx
VHAnalyser::VHAnalyser(TTree *tree, std::string fname)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MC_RelValTTbar.root");
      if (!f) {
         f = new TFile("MC_RelValTTbar.root");
         f->cd("MC_RelValTTbar.root:/treeCreator");
      }
      tree = (TTree*)gDirectory->Get("vhtree");

   }
   Init(tree);
   fileName = fname;
}

VHAnalyser::~VHAnalyser()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t VHAnalyser::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t VHAnalyser::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void VHAnalyser::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_, &b_Event_);
   fChain->SetBranchAddress("Event.fUniqueID", Event_fUniqueID, &b_Event_fUniqueID);
   fChain->SetBranchAddress("Event.fBits", Event_fBits, &b_Event_fBits);
   fChain->SetBranchAddress("Event.run", Event_run, &b_Event_run);
   fChain->SetBranchAddress("Event.event", Event_event, &b_Event_event);
   fChain->SetBranchAddress("Event.lumis", Event_lumis, &b_Event_lumis);
   fChain->SetBranchAddress("Event.bunch", Event_bunch, &b_Event_bunch);
   fChain->SetBranchAddress("Event.orbit", Event_orbit, &b_Event_orbit);
   fChain->SetBranchAddress("Event.time", Event_time, &b_Event_time);
   fChain->SetBranchAddress("Event.isdata", Event_isdata, &b_Event_isdata);
   fChain->SetBranchAddress("Vertex", &Vertex_, &b_Vertex_);
   fChain->SetBranchAddress("Vertex.fUniqueID", Vertex_fUniqueID, &b_Vertex_fUniqueID);
   fChain->SetBranchAddress("Vertex.fBits", Vertex_fBits, &b_Vertex_fBits);
   fChain->SetBranchAddress("Vertex.x", Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex.y", Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex.z", Vertex_z, &b_Vertex_z);
   fChain->SetBranchAddress("Vertex.xErr", Vertex_xErr, &b_Vertex_xErr);
   fChain->SetBranchAddress("Vertex.yErr", Vertex_yErr, &b_Vertex_yErr);
   fChain->SetBranchAddress("Vertex.zErr", Vertex_zErr, &b_Vertex_zErr);
   fChain->SetBranchAddress("Vertex.rho", Vertex_rho, &b_Vertex_rho);
   fChain->SetBranchAddress("Vertex.chi2", Vertex_chi2, &b_Vertex_chi2);
   fChain->SetBranchAddress("Vertex.ndf", Vertex_ndf, &b_Vertex_ndf);
   fChain->SetBranchAddress("Vertex.ntracks", Vertex_ntracks, &b_Vertex_ntracks);
   fChain->SetBranchAddress("Vertex.ntracksw05", Vertex_ntracksw05, &b_Vertex_ntracksw05);
   fChain->SetBranchAddress("Vertex.isfake", Vertex_isfake, &b_Vertex_isfake);
   fChain->SetBranchAddress("nVertex", &nVertex, &b_fnVertex);
   fChain->SetBranchAddress("CaloJet", &CaloJet_, &b_CaloJet_);
   fChain->SetBranchAddress("CaloJet.fUniqueID", CaloJet_fUniqueID, &b_CaloJet_fUniqueID);
   fChain->SetBranchAddress("CaloJet.fBits", CaloJet_fBits, &b_CaloJet_fBits);
   fChain->SetBranchAddress("CaloJet.eta", CaloJet_eta, &b_CaloJet_eta);
   fChain->SetBranchAddress("CaloJet.phi", CaloJet_phi, &b_CaloJet_phi);
   fChain->SetBranchAddress("CaloJet.pt", CaloJet_pt, &b_CaloJet_pt);
   fChain->SetBranchAddress("CaloJet.pt_raw", CaloJet_pt_raw, &b_CaloJet_pt_raw);
   fChain->SetBranchAddress("CaloJet.energy", CaloJet_energy, &b_CaloJet_energy);
   fChain->SetBranchAddress("CaloJet.energy_raw", CaloJet_energy_raw, &b_CaloJet_energy_raw);
   fChain->SetBranchAddress("CaloJet.jecUnc", CaloJet_jecUnc, &b_CaloJet_jecUnc);
   fChain->SetBranchAddress("CaloJet.resJEC", CaloJet_resJEC, &b_CaloJet_resJEC);
   fChain->SetBranchAddress("CaloJet.overlaps", CaloJet_overlaps, &b_CaloJet_overlaps);
   fChain->SetBranchAddress("CaloJet.partonFlavour", CaloJet_partonFlavour, &b_CaloJet_partonFlavour);
   fChain->SetBranchAddress("CaloJet.emf", CaloJet_emf, &b_CaloJet_emf);
   fChain->SetBranchAddress("CaloJet.resEmf", CaloJet_resEmf, &b_CaloJet_resEmf);
   fChain->SetBranchAddress("CaloJet.hadf", CaloJet_hadf, &b_CaloJet_hadf);
   fChain->SetBranchAddress("CaloJet.n90Hits", CaloJet_n90Hits, &b_CaloJet_n90Hits);
   fChain->SetBranchAddress("CaloJet.fHPD", CaloJet_fHPD, &b_CaloJet_fHPD);
   fChain->SetBranchAddress("CaloJet.fRBX", CaloJet_fRBX, &b_CaloJet_fRBX);
   fChain->SetBranchAddress("CaloJet.sigmaEta", CaloJet_sigmaEta, &b_CaloJet_sigmaEta);
   fChain->SetBranchAddress("CaloJet.sigmaPhi", CaloJet_sigmaPhi, &b_CaloJet_sigmaPhi);
   fChain->SetBranchAddress("CaloJet.trackCountingHighEffBTag", CaloJet_trackCountingHighEffBTag, &b_CaloJet_trackCountingHighEffBTag);
   fChain->SetBranchAddress("CaloJet.trackCountingHighPurBTag", CaloJet_trackCountingHighPurBTag, &b_CaloJet_trackCountingHighPurBTag);
   fChain->SetBranchAddress("CaloJet.simpleSecondaryVertexHighEffBTag", CaloJet_simpleSecondaryVertexHighEffBTag, &b_CaloJet_simpleSecondaryVertexHighEffBTag);
   fChain->SetBranchAddress("CaloJet.simpleSecondaryVertexHighPurBTag", CaloJet_simpleSecondaryVertexHighPurBTag, &b_CaloJet_simpleSecondaryVertexHighPurBTag);
   fChain->SetBranchAddress("CaloJet.jetProbabilityBTag", CaloJet_jetProbabilityBTag, &b_CaloJet_jetProbabilityBTag);
   fChain->SetBranchAddress("CaloJet.jetBProbabilityBTag", CaloJet_jetBProbabilityBTag, &b_CaloJet_jetBProbabilityBTag);
   fChain->SetBranchAddress("CaloJet.passLooseID", CaloJet_passLooseID, &b_CaloJet_passLooseID);
   fChain->SetBranchAddress("CaloJet.passTightID", CaloJet_passTightID, &b_CaloJet_passTightID);
   fChain->SetBranchAddress("nCaloJet", &nCaloJet, &b_fnCaloJet);
   fChain->SetBranchAddress("Jet", &Jet_, &b_Jet_);
   fChain->SetBranchAddress("Jet.fUniqueID", &Jet_fUniqueID, &b_Jet_fUniqueID);
   fChain->SetBranchAddress("Jet.fBits", &Jet_fBits, &b_Jet_fBits);
   fChain->SetBranchAddress("Jet.eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet.phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet.pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet.pt_raw", &Jet_pt_raw, &b_Jet_pt_raw);
   fChain->SetBranchAddress("Jet.energy", &Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet.energy_raw", &Jet_energy_raw, &b_Jet_energy_raw);
   fChain->SetBranchAddress("Jet.jecUnc", &Jet_jecUnc, &b_Jet_jecUnc);
   fChain->SetBranchAddress("Jet.resJEC", &Jet_resJEC, &b_Jet_resJEC);
   fChain->SetBranchAddress("Jet.partonFlavour", &Jet_partonFlavour, &b_Jet_partonFlavour);
   fChain->SetBranchAddress("Jet.chargedEmEnergyFraction", &Jet_chargedEmEnergyFraction, &b_Jet_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jet.chargedHadronEnergyFraction", &Jet_chargedHadronEnergyFraction, &b_Jet_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jet.chargedMuEnergyFraction", &Jet_chargedMuEnergyFraction, &b_Jet_chargedMuEnergyFraction);
   fChain->SetBranchAddress("Jet.electronEnergyFraction", &Jet_electronEnergyFraction, &b_Jet_electronEnergyFraction);
   fChain->SetBranchAddress("Jet.muonEnergyFraction", &Jet_muonEnergyFraction, &b_Jet_muonEnergyFraction);
   fChain->SetBranchAddress("Jet.neutralEmEnergyFraction", &Jet_neutralEmEnergyFraction, &b_Jet_neutralEmEnergyFraction);
   fChain->SetBranchAddress("Jet.neutralHadronEnergyFraction", &Jet_neutralHadronEnergyFraction, &b_Jet_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("Jet.photonEnergyFraction", &Jet_photonEnergyFraction, &b_Jet_photonEnergyFraction);
   fChain->SetBranchAddress("Jet.chargedHadronMultiplicity", &Jet_chargedHadronMultiplicity, &b_Jet_chargedHadronMultiplicity);
   fChain->SetBranchAddress("Jet.chargedMultiplicity", &Jet_chargedMultiplicity, &b_Jet_chargedMultiplicity);
   fChain->SetBranchAddress("Jet.electronMultiplicity", &Jet_electronMultiplicity, &b_Jet_electronMultiplicity);
   fChain->SetBranchAddress("Jet.muonMultiplicity", &Jet_muonMultiplicity, &b_Jet_muonMultiplicity);
   fChain->SetBranchAddress("Jet.neutralHadronMultiplicity", &Jet_neutralHadronMultiplicity, &b_Jet_neutralHadronMultiplicity);
   fChain->SetBranchAddress("Jet.neutralMultiplicity", &Jet_neutralMultiplicity, &b_Jet_neutralMultiplicity);
   fChain->SetBranchAddress("Jet.photonMultiplicity", &Jet_photonMultiplicity, &b_Jet_photonMultiplicity);
   fChain->SetBranchAddress("Jet.nConstituents", &Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet.trackCountingHighEffBTag", &Jet_trackCountingHighEffBTag, &b_Jet_trackCountingHighEffBTag);
   fChain->SetBranchAddress("Jet.trackCountingHighPurBTag", &Jet_trackCountingHighPurBTag, &b_Jet_trackCountingHighPurBTag);
   fChain->SetBranchAddress("Jet.simpleSecondaryVertexHighEffBTag", &Jet_simpleSecondaryVertexHighEffBTag, &b_Jet_simpleSecondaryVertexHighEffBTag);
   fChain->SetBranchAddress("Jet.simpleSecondaryVertexHighPurBTag", &Jet_simpleSecondaryVertexHighPurBTag, &b_Jet_simpleSecondaryVertexHighPurBTag);
   fChain->SetBranchAddress("Jet.jetProbabilityBTag", &Jet_jetProbabilityBTag, &b_Jet_jetProbabilityBTag);
   fChain->SetBranchAddress("Jet.jetBProbabilityBTag", &Jet_jetBProbabilityBTag, &b_Jet_jetBProbabilityBTag);
   fChain->SetBranchAddress("Jet.passLooseID", &Jet_passLooseID, &b_Jet_passLooseID);
   fChain->SetBranchAddress("Jet.passTightID", &Jet_passTightID, &b_Jet_passTightID);
   fChain->SetBranchAddress("nJet", &nJet, &b_fnJet);
   fChain->SetBranchAddress("Electron", &Electron_, &b_Electron_);
   fChain->SetBranchAddress("Electron.fUniqueID", Electron_fUniqueID, &b_Electron_fUniqueID);
   fChain->SetBranchAddress("Electron.fBits", Electron_fBits, &b_Electron_fBits);
   fChain->SetBranchAddress("Electron.eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron.phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron.pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron.trackPt", Electron_trackPt, &b_Electron_trackPt);
   fChain->SetBranchAddress("Electron.energy", Electron_energy, &b_Electron_energy);
   fChain->SetBranchAddress("Electron.caloEnergy", Electron_caloEnergy, &b_Electron_caloEnergy);
   fChain->SetBranchAddress("Electron.charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron.hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron.sigmaEtaEta", Electron_sigmaEtaEta, &b_Electron_sigmaEtaEta);
   fChain->SetBranchAddress("Electron.sigmaIEtaIEta", Electron_sigmaIEtaIEta, &b_Electron_sigmaIEtaIEta);
   fChain->SetBranchAddress("Electron.deltaPhiTrkSC", Electron_deltaPhiTrkSC, &b_Electron_deltaPhiTrkSC);
   fChain->SetBranchAddress("Electron.deltaEtaTrkSC", Electron_deltaEtaTrkSC, &b_Electron_deltaEtaTrkSC);
   fChain->SetBranchAddress("Electron.classif", Electron_classif, &b_Electron_classif);
   fChain->SetBranchAddress("Electron.e1x5overe5x5", Electron_e1x5overe5x5, &b_Electron_e1x5overe5x5);
   fChain->SetBranchAddress("Electron.e2x5overe5x5", Electron_e2x5overe5x5, &b_Electron_e2x5overe5x5);
   fChain->SetBranchAddress("Electron.isoEcal03", Electron_isoEcal03, &b_Electron_isoEcal03);
   fChain->SetBranchAddress("Electron.isoHcal03", Electron_isoHcal03, &b_Electron_isoHcal03);
   fChain->SetBranchAddress("Electron.isoTrk03", Electron_isoTrk03, &b_Electron_isoTrk03);
   fChain->SetBranchAddress("Electron.isoEcal04", Electron_isoEcal04, &b_Electron_isoEcal04);
   fChain->SetBranchAddress("Electron.isoHcal04", Electron_isoHcal04, &b_Electron_isoHcal04);
   fChain->SetBranchAddress("Electron.isoTrk04", Electron_isoTrk04, &b_Electron_isoTrk04);
   fChain->SetBranchAddress("Electron.isoRel03", Electron_isoRel03, &b_Electron_isoRel03);
   fChain->SetBranchAddress("Electron.isoRel04", Electron_isoRel04, &b_Electron_isoRel04);
   fChain->SetBranchAddress("Electron.missingHits", Electron_missingHits, &b_Electron_missingHits);
   fChain->SetBranchAddress("Electron.dist_vec", Electron_dist_vec, &b_Electron_dist_vec);
   fChain->SetBranchAddress("Electron.dCotTheta", Electron_dCotTheta, &b_Electron_dCotTheta);
   fChain->SetBranchAddress("Electron.scEn", Electron_scEn, &b_Electron_scEn);
   fChain->SetBranchAddress("Electron.scEta", Electron_scEta, &b_Electron_scEta);
   fChain->SetBranchAddress("Electron.scPhi", Electron_scPhi, &b_Electron_scPhi);
   fChain->SetBranchAddress("Electron.scET", Electron_scET, &b_Electron_scET);
   fChain->SetBranchAddress("Electron.scRawEnergy", Electron_scRawEnergy, &b_Electron_scRawEnergy);
   fChain->SetBranchAddress("Electron.vtxDist3D", Electron_vtxDist3D, &b_Electron_vtxDist3D);
   fChain->SetBranchAddress("Electron.vtxIndex", Electron_vtxIndex, &b_Electron_vtxIndex);
   fChain->SetBranchAddress("Electron.vtxDistZ", Electron_vtxDistZ, &b_Electron_vtxDistZ);
   fChain->SetBranchAddress("Electron.pfRelIso", Electron_pfRelIso, &b_Electron_pfRelIso);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_fnElectron);
   fChain->SetBranchAddress("MET", &MET_, &b_MET_);
   fChain->SetBranchAddress("MET.fUniqueID", MET_fUniqueID, &b_MET_fUniqueID);
   fChain->SetBranchAddress("MET.fBits", MET_fBits, &b_MET_fBits);
   fChain->SetBranchAddress("MET.met", MET_met, &b_MET_met);
   fChain->SetBranchAddress("MET.metphi", MET_metphi, &b_MET_metphi);
   fChain->SetBranchAddress("MET.sumet", MET_sumet, &b_MET_sumet);
   fChain->SetBranchAddress("MET.metuncorr", MET_metuncorr, &b_MET_metuncorr);
   fChain->SetBranchAddress("MET.metphiuncorr", MET_metphiuncorr, &b_MET_metphiuncorr);
   fChain->SetBranchAddress("MET.sumetuncorr", MET_sumetuncorr, &b_MET_sumetuncorr);
   fChain->SetBranchAddress("nMET", &nMET, &b_fnMET);
   fChain->SetBranchAddress("Muon", &Muon_, &b_Muon_);
   fChain->SetBranchAddress("Muon.fUniqueID", Muon_fUniqueID, &b_Muon_fUniqueID);
   fChain->SetBranchAddress("Muon.fBits", Muon_fBits, &b_Muon_fBits);
   fChain->SetBranchAddress("Muon.eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon.phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon.pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon.p", Muon_p, &b_Muon_p);
   fChain->SetBranchAddress("Muon.energy", Muon_energy, &b_Muon_energy);
   fChain->SetBranchAddress("Muon.charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon.trkD0", Muon_trkD0, &b_Muon_trkD0);
   fChain->SetBranchAddress("Muon.trkD0Error", Muon_trkD0Error, &b_Muon_trkD0Error);
   fChain->SetBranchAddress("Muon.trkDz", Muon_trkDz, &b_Muon_trkDz);
   fChain->SetBranchAddress("Muon.trkDzError", Muon_trkDzError, &b_Muon_trkDzError);
   fChain->SetBranchAddress("Muon.globalChi2", Muon_globalChi2, &b_Muon_globalChi2);
   fChain->SetBranchAddress("Muon.trkIso", Muon_trkIso, &b_Muon_trkIso);
   fChain->SetBranchAddress("Muon.ecalIso", Muon_ecalIso, &b_Muon_ecalIso);
   fChain->SetBranchAddress("Muon.hcalIso", Muon_hcalIso, &b_Muon_hcalIso);
   fChain->SetBranchAddress("Muon.hoIso", Muon_hoIso, &b_Muon_hoIso);
   fChain->SetBranchAddress("Muon.relIso", Muon_relIso, &b_Muon_relIso);
   fChain->SetBranchAddress("Muon.passID", Muon_passID, &b_Muon_passID);
   fChain->SetBranchAddress("Muon.vtxDist3D", Muon_vtxDist3D, &b_Muon_vtxDist3D);
   fChain->SetBranchAddress("Muon.vtxIndex", Muon_vtxIndex, &b_Muon_vtxIndex);
   fChain->SetBranchAddress("Muon.vtxDistZ", Muon_vtxDistZ, &b_Muon_vtxDistZ);
   fChain->SetBranchAddress("Muon.pixHits", Muon_pixHits, &b_Muon_pixHits);
   fChain->SetBranchAddress("Muon.trkHits", Muon_trkHits, &b_Muon_trkHits);
   fChain->SetBranchAddress("Muon.matches", Muon_matches, &b_Muon_matches);
   fChain->SetBranchAddress("Muon.pfRelIso", Muon_pfRelIso, &b_Muon_pfRelIso);
   fChain->SetBranchAddress("Muon.isTrackerMuon", Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_fnMuon);
   fChain->SetBranchAddress("Tau", &Tau_, &b_Tau_);
   fChain->SetBranchAddress("Tau.fUniqueID", Tau_fUniqueID, &b_Tau_fUniqueID);
   fChain->SetBranchAddress("Tau.fBits", Tau_fBits, &b_Tau_fBits);
   fChain->SetBranchAddress("Tau.eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau.phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau.pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau.energy", Tau_energy, &b_Tau_energy);
   fChain->SetBranchAddress("Tau.charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau.mass", Tau_mass, &b_Tau_mass);
   fChain->SetBranchAddress("Tau.leadChargedParticlePt", Tau_leadChargedParticlePt, &b_Tau_leadChargedParticlePt);
   fChain->SetBranchAddress("Tau.leadNeutralParticlePt", Tau_leadNeutralParticlePt, &b_Tau_leadNeutralParticlePt);
   fChain->SetBranchAddress("Tau.leadParticlePt", Tau_leadParticlePt, &b_Tau_leadParticlePt);
   fChain->SetBranchAddress("Tau.numChargedHadronsSignalCone", Tau_numChargedHadronsSignalCone, &b_Tau_numChargedHadronsSignalCone);
   fChain->SetBranchAddress("Tau.numNeutralHadronsSignalCone", Tau_numNeutralHadronsSignalCone, &b_Tau_numNeutralHadronsSignalCone);
   fChain->SetBranchAddress("Tau.numPhotonsSignalCone", Tau_numPhotonsSignalCone, &b_Tau_numPhotonsSignalCone);
   fChain->SetBranchAddress("Tau.numParticlesSignalCone", Tau_numParticlesSignalCone, &b_Tau_numParticlesSignalCone);
   fChain->SetBranchAddress("Tau.numChargedHadronsIsoCone", Tau_numChargedHadronsIsoCone, &b_Tau_numChargedHadronsIsoCone);
   fChain->SetBranchAddress("Tau.numNeutralHadronsIsoCone", Tau_numNeutralHadronsIsoCone, &b_Tau_numNeutralHadronsIsoCone);
   fChain->SetBranchAddress("Tau.numPhotonsIsoCone", Tau_numPhotonsIsoCone, &b_Tau_numPhotonsIsoCone);
   fChain->SetBranchAddress("Tau.numParticlesIsoCone", Tau_numParticlesIsoCone, &b_Tau_numParticlesIsoCone);
   fChain->SetBranchAddress("Tau.ptSumPFChargedHadronsIsoCone", Tau_ptSumPFChargedHadronsIsoCone, &b_Tau_ptSumPFChargedHadronsIsoCone);
   fChain->SetBranchAddress("Tau.ptSumPhotonsIsoCone", Tau_ptSumPhotonsIsoCone, &b_Tau_ptSumPhotonsIsoCone);
   fChain->SetBranchAddress("Tau.decayModeFinding", Tau_decayModeFinding, &b_Tau_decayModeFinding);
   fChain->SetBranchAddress("Tau.looseIsolation", Tau_looseIsolation, &b_Tau_looseIsolation);
   fChain->SetBranchAddress("Tau.mediumIsolation", Tau_mediumIsolation, &b_Tau_mediumIsolation);
   fChain->SetBranchAddress("Tau.tightIsolation", Tau_tightIsolation, &b_Tau_tightIsolation);
   fChain->SetBranchAddress("Tau.againstMuonLoose", Tau_againstMuonLoose, &b_Tau_againstMuonLoose);
   fChain->SetBranchAddress("Tau.againstMuonTight", Tau_againstMuonTight, &b_Tau_againstMuonTight);
   fChain->SetBranchAddress("Tau.againstElectronLoose", Tau_againstElectronLoose, &b_Tau_againstElectronLoose);
   fChain->SetBranchAddress("Tau.againstElectronMedium", Tau_againstElectronMedium, &b_Tau_againstElectronMedium);
   fChain->SetBranchAddress("Tau.againstElectronTight", Tau_againstElectronTight, &b_Tau_againstElectronTight);
   fChain->SetBranchAddress("Tau.pfElectronMVA", Tau_pfElectronMVA, &b_Tau_pfElectronMVA);
   fChain->SetBranchAddress("Tau.jetPt", Tau_jetPt, &b_Tau_jetPt);
   fChain->SetBranchAddress("Tau.jetEta", Tau_jetEta, &b_Tau_jetEta);
   fChain->SetBranchAddress("Tau.jetPhi", Tau_jetPhi, &b_Tau_jetPhi);
   fChain->SetBranchAddress("nTau", &nTau, &b_fnTau);
   fChain->SetBranchAddress("GenParticle", &GenParticle_, &b_GenParticle_);
   fChain->SetBranchAddress("GenParticle.fUniqueID", GenParticle_fUniqueID, &b_GenParticle_fUniqueID);
   fChain->SetBranchAddress("GenParticle.fBits", GenParticle_fBits, &b_GenParticle_fBits);
   fChain->SetBranchAddress("GenParticle.eta", GenParticle_eta, &b_GenParticle_eta);
   fChain->SetBranchAddress("GenParticle.phi", GenParticle_phi, &b_GenParticle_phi);
   fChain->SetBranchAddress("GenParticle.p", GenParticle_p, &b_GenParticle_p);
   fChain->SetBranchAddress("GenParticle.px", GenParticle_px, &b_GenParticle_px);
   fChain->SetBranchAddress("GenParticle.py", GenParticle_py, &b_GenParticle_py);
   fChain->SetBranchAddress("GenParticle.pz", GenParticle_pz, &b_GenParticle_pz);
   fChain->SetBranchAddress("GenParticle.pt", GenParticle_pt, &b_GenParticle_pt);
   fChain->SetBranchAddress("GenParticle.energy", GenParticle_energy, &b_GenParticle_energy);
   fChain->SetBranchAddress("GenParticle.pdgId", GenParticle_pdgId, &b_GenParticle_pdgId);
   fChain->SetBranchAddress("GenParticle.vx", GenParticle_vx, &b_GenParticle_vx);
   fChain->SetBranchAddress("GenParticle.vy", GenParticle_vy, &b_GenParticle_vy);
   fChain->SetBranchAddress("GenParticle.vz", GenParticle_vz, &b_GenParticle_vz);
   fChain->SetBranchAddress("GenParticle.numDaught", GenParticle_numDaught, &b_GenParticle_numDaught);
   fChain->SetBranchAddress("GenParticle.status", GenParticle_status, &b_GenParticle_status);
   fChain->SetBranchAddress("GenParticle.motherIndex", GenParticle_motherIndex, &b_GenParticle_motherIndex);
   fChain->SetBranchAddress("nGenParticle", &nGenParticle, &b_fnGenParticle);
   fChain->SetBranchAddress("GenJet", &GenJet_, &b_GenJet_);
   fChain->SetBranchAddress("GenJet.fUniqueID", GenJet_fUniqueID, &b_GenJet_fUniqueID);
   fChain->SetBranchAddress("GenJet.fBits", GenJet_fBits, &b_GenJet_fBits);
   fChain->SetBranchAddress("GenJet.eta", GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet.phi", GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet.p", GenJet_p, &b_GenJet_p);
   fChain->SetBranchAddress("GenJet.pt", GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("GenJet.energy", GenJet_energy, &b_GenJet_energy);
   fChain->SetBranchAddress("GenJet.emf", GenJet_emf, &b_GenJet_emf);
   fChain->SetBranchAddress("GenJet.hadf", GenJet_hadf, &b_GenJet_hadf);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_fnGenJet);
   fChain->SetBranchAddress("GenMET", &GenMET_, &b_GenMET_);
   fChain->SetBranchAddress("GenMET.fUniqueID", GenMET_fUniqueID, &b_GenMET_fUniqueID);
   fChain->SetBranchAddress("GenMET.fBits", GenMET_fBits, &b_GenMET_fBits);
   fChain->SetBranchAddress("GenMET.met", GenMET_met, &b_GenMET_met);
   fChain->SetBranchAddress("GenMET.metphi", GenMET_metphi, &b_GenMET_metphi);
   fChain->SetBranchAddress("GenMET.sumet", GenMET_sumet, &b_GenMET_sumet);
   fChain->SetBranchAddress("nGenMET", &nGenMET, &b_fnGenMET);
   Notify();
}

Bool_t VHAnalyser::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void VHAnalyser::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t VHAnalyser::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef VHAnalyser_cxx
