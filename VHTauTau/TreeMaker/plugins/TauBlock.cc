#include <iostream>
#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Utilities/General/interface/FileInPath.h"
#include "Bianchi/Utilities/interface/AntiElectronIDMVA.h"

#include "VHTauTau/TreeMaker/plugins/TauBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

TauBlock::TauBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("patTauSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc"))
{
  std::string method = iConfig.getParameter<std::string>("methodName");
  edm::FileInPath weightsX0BL = iConfig.getParameter<edm::FileInPath>("weightsX0BL");
  edm::FileInPath weights11BL = iConfig.getParameter<edm::FileInPath>("weights11BL");
  edm::FileInPath weights01BL = iConfig.getParameter<edm::FileInPath>("weights01BL");
  edm::FileInPath weightsX0EC = iConfig.getParameter<edm::FileInPath>("weightsX0EC");
  edm::FileInPath weights11EC = iConfig.getParameter<edm::FileInPath>("weights11EC");
  edm::FileInPath weights01EC = iConfig.getParameter<edm::FileInPath>("weights01EC");

  antiE = new AntiElectronIDMVA();
  antiE->Initialize(method,
                    weightsX0BL.fullPath(), 
                    weights11BL.fullPath(), 
                    weights01BL.fullPath(), 
                    weightsX0EC.fullPath(), 
                    weights11EC.fullPath(), 
                    weights01EC.fullPath() 
  );
}
TauBlock::~TauBlock() { delete antiE; }
void TauBlock::beginJob() 
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneTau = new TClonesArray("vhtm::Tau");
  tree->Branch("Tau", &cloneTau, 32000, 2);
  tree->Branch("nTau", &fnTau, "fnTau/I");
}
void TauBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneTau->Clear();
  fnTau = 0;

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(_inputTag, taus);
  
  if (taus.isValid()) {
    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    
    for (std::vector<pat::Tau>::const_iterator it  = taus->begin(); 
                                               it != taus->end(); ++it) {
      if (fnTau == kMaxTau) {
	edm::LogInfo("TauBlock") << "Too many PAT Taus, fnTau = " << fnTau;
	break;
      }
      tauB = new ((*cloneTau)[fnTau++]) vhtm::Tau();

      // Store Tau variables
      tauB->eta    = it->eta();
      tauB->phi    = it->phi();
      tauB->pt     = it->pt();
      tauB->energy = it->energy();
      tauB->charge = it->charge();
      if (it->leadTrack().isAvailable() && it->leadTrack().isNonnull()) {
	reco::TrackRef trk = it->leadTrack();
        tauB->leadTrkPt     = trk->pt();
        tauB->leadTrkEta    = trk->eta();
        tauB->leadTrkPhi    = trk->phi();
        tauB->leadTrkCharge = trk->charge();
      }

      // IP of leadPFChargedHadrCand wrt PV
      // Vertex association
      double minVtxDist3D = 9999.;
      int indexVtx = -1;
      double vertexDz = 9999.;
      double vertexDxy = 9999.;
      if (primaryVertices.isValid() && it->leadTrack().isAvailable() && it->leadTrack().isNonnull()) {
	edm::LogInfo("TauBlock") << "Total # Primary Vertices: " << primaryVertices->size();

        for (reco::VertexCollection::const_iterator v_it  = primaryVertices->begin(); 
                                                    v_it != primaryVertices->end(); ++v_it) {
          double dxy = it->leadTrack()->dxy(v_it->position());
          double dz  = it->leadTrack()->dz(v_it->position());
          double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
          if (dist3D < minVtxDist3D) {
            minVtxDist3D = dist3D;
            indexVtx = int(std::distance(primaryVertices->begin(), v_it));
            vertexDxy = dxy;
            vertexDz = dz;
          }
        }
      } 
      else {
	edm::LogError("TauBlock") << "Error >> Failed to get VertexCollection for label: " 
                                  << _vtxInputTag;
      }
      tauB->vtxIndex = indexVtx;
      tauB->vtxDxy   = vertexDxy;
      tauB->vtxDz    = vertexDz;

      // Leading particle pT
      tauB->leadChargedParticlePt = it->leadPFChargedHadrCand().isNonnull() 
                                      ? it->leadPFChargedHadrCand()->pt(): 0.;
      tauB->leadNeutralParticlePt = it->leadPFNeutralCand().isNonnull() 
                                      ? it->leadPFNeutralCand()->et(): 0.;
      tauB->leadParticlePt        = it->leadPFCand().isNonnull() 
                                      ? it->leadPFCand()->et(): 0.;      
      // Number of charged/neutral candidates and photons in different cones
      tauB->numChargedHadronsSignalCone = it->signalPFChargedHadrCands().size();
      tauB->numNeutralHadronsSignalCone = it->signalPFNeutrHadrCands().size();
      tauB->numPhotonsSignalCone        = it->signalPFGammaCands().size();
      tauB->numParticlesSignalCone      = it->signalPFCands().size();
      //tauB->numPi0SignalCone            = it->signalPiZeroCandidates().size();
      
      tauB->numChargedHadronsIsoCone = it->isolationPFChargedHadrCands().size();
      tauB->numNeutralHadronsIsoCone = it->isolationPFNeutrHadrCands().size();
      tauB->numPhotonsIsoCone        = it->isolationPFGammaCands().size();
      tauB->numParticlesIsoCone      = it->isolationPFCands().size();
      //tauB->numPi0IsoCone            = it->isolationPiZeroCandidates().size();
      
      tauB->ptSumPFChargedHadronsIsoCone = it->isolationPFChargedHadrCandsPtSum();
      tauB->ptSumPhotonsIsoCone          = it->isolationPFGammaCandsEtSum();
      const reco::PFCandidateRefVector& nvList = it->isolationPFNeutrHadrCands();
      double ptSum = 0;
      for (reco::PFCandidateRefVector::const_iterator iCand  = nvList.begin(); 
                                                      iCand != nvList.end(); ++iCand) {
	ptSum += std::abs((**iCand).pt());
      }  
      tauB->ptSumPFNeutralHadronsIsoCone = ptSum;

      // tau id. discriminators
      tauB->decayModeFinding = it->tauID("decayModeFinding");
      tauB->looseIsolation   = it->tauID("byLooseIsolation");
      tauB->mediumIsolation  = it->tauID("byMediumIsolation");
      tauB->tightIsolation   = it->tauID("byTightIsolation");

      // discriminators against electrons/muons
      tauB->againstMuonLoose      = it->tauID("againstMuonLoose");
      tauB->againstMuonTight      = it->tauID("againstMuonTight");
      tauB->againstElectronLoose  = it->tauID("againstElectronLoose");
      tauB->againstElectronMedium = it->tauID("againstElectronMedium");
      tauB->againstElectronTight  = it->tauID("againstElectronTight");
      tauB->pfElectronMVA         = it->leadPFCand().isNonnull() 
                                  ? it->leadPFCand()->mva_e_pi() : 1.;

      tauB->byVLooseCombinedIsolationDeltaBetaCorr = it->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
      tauB->byLooseCombinedIsolationDeltaBetaCorr  = it->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      tauB->byMediumCombinedIsolationDeltaBetaCorr = it->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      tauB->byTightCombinedIsolationDeltaBetaCorr  = it->tauID("byTightCombinedIsolationDeltaBetaCorr");
      tauB->byVLooseIsolationDeltaBetaCorr         = it->tauID("byVLooseIsolationDeltaBetaCorr");
      tauB->byLooseIsolationDeltaBetaCorr          = it->tauID("byLooseIsolationDeltaBetaCorr");
      tauB->byMediumIsolationDeltaBetaCorr         = it->tauID("byMediumIsolationDeltaBetaCorr");
      tauB->byTightIsolationDeltaBetaCorr          = it->tauID("byTightIsolationDeltaBetaCorr");

      // kinematic variables for PFJet associated to PFTau
      if (it->pfJetRef().isAvailable() && it->pfJetRef().isNonnull()) {
	const reco::PFJetRef &jtr = it->pfJetRef();
        tauB->jetPt  = jtr->pt();
        tauB->jetEta = jtr->eta();
        tauB->jetPhi = jtr->phi();
      }

      // NEW quantities
      tauB->emFraction              = it->emFraction(); 
      tauB->maximumHCALPFClusterEt  = it->maximumHCALPFClusterEt();
      tauB->ecalStripSumEOverPLead  = it->ecalStripSumEOverPLead();
      tauB->bremsRecoveryEOverPLead = it->bremsRecoveryEOverPLead();
      tauB->hcalTotOverPLead        = it->hcalTotOverPLead();
      tauB->hcalMaxOverPLead        = it->hcalMaxOverPLead();
      tauB->hcal3x3OverPLead        = it->hcal3x3OverPLead();

      tauB->etaetaMoment = it->etaetaMoment();
      tauB->phiphiMoment = it->phiphiMoment();
      tauB->phiphiMoment = it->etaphiMoment();
      
      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      tauB->vx = vertex.x();             
      tauB->vy = vertex.y();             
      tauB->vz = vertex.z();             

      tauB->zvertex = it->vz(); // distance from the primary vertex
      tauB->mass    = it->p4().M();
      tauB->ltsipt  = TMath::Abs(it->leadPFChargedHadrCandsignedSipt());

      // ElectronIDMVA, electron faking tau
      const pat::Tau& tau = *it;  
      tauB->mva = (it->leadPFChargedHadrCand().isNonnull()) ? antiE->MVAValue(&tau) : -1;
      tauB->againstElectronMVA = it->tauID("againstElectronMVA");
    }
  }
  else {
    edm::LogError("TauBlock") << "Error! Failed to get pat::Tau collection for label: " 
                              << _inputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
