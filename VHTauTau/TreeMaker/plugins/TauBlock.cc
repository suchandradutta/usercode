#include "TTree.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/TauBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

TauBlock::TauBlock(const edm::ParameterSet& iConfig) :
  _tree(0),
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("patTauSrc"))
{}
void TauBlock::beginJob() 
{
  // Get TTree pointer
  //edm::Service<TFileService> fs;
  //TTree* tree = fs->getObject<TTree>("vhtree");
  if (!_tree) _tree = Utility::getTree("vhtree");
  cloneTau = new TClonesArray("Tau");
  _tree->Branch("Tau", &cloneTau, 32000, 2);
  _tree->Branch("nTau", &fnTau, "fnTau/I");
}
void TauBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneTau->Clear();
  fnTau = 0;

  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(_inputTag, taus);
  
  if (taus.isValid()) {
    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    
    for (std::vector<pat::Tau>::const_iterator it = taus->begin(); 
                                              it != taus->end(); ++it ) {
      
      if (fnTau == kMaxTau) {
	edm::LogInfo("TauBlock") << "Too many PAT Taus, fnTau = " << fnTau;
	break;
      }
      tauB = new ((*cloneTau)[fnTau++]) Tau();
      // Store Tau variables
      tauB->eta    = it->eta();
      tauB->phi    = it->phi();
      tauB->pt     = it->pt();
      tauB->energy = it->energy();
      tauB->charge = it->charge();
      tauB->mass   = it->mass();

      // Leading particle pT
      tauB->leadChargedParticlePt = it->leadPFChargedHadrCand().isNonnull() ? it->leadPFChargedHadrCand()->et(): 0.;
      tauB->leadNeutralParticlePt = it->leadPFNeutralCand().isNonnull() ? it->leadPFNeutralCand()->et(): 0.;
      tauB->leadParticlePt        = it->leadPFCand().isNonnull() ? it->leadPFCand()->et(): 0.;      

      // Number of charged/neutral candidates and photons in different cones
      tauB->numChargedHadronsSignalCone = it->signalPFChargedHadrCands().size();
      tauB->numNeutralHadronsSignalCone = it->signalPFNeutrHadrCands().size();
      tauB->numPhotonsSignalCone        = it->signalPFGammaCands().size();
      tauB->numParticlesSignalCone      = it->signalPFCands().size();
      
      tauB->numChargedHadronsIsoCone = it->isolationPFChargedHadrCands().size();
      tauB->numNeutralHadronsIsoCone = it->isolationPFNeutrHadrCands().size();
      tauB->numPhotonsIsoCone        = it->isolationPFGammaCands().size();
      tauB->numParticlesIsoCone      = it->isolationPFCands().size();
      
      tauB->ptSumPFChargedHadronsIsoCone = it->isolationPFChargedHadrCandsPtSum();
      tauB->ptSumPhotonsIsoCone          = it->isolationPFGammaCandsEtSum();

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
      tauB->pfElectronMVA         = it->leadPFCand().isNonnull() ? it->leadPFCand()->mva_e_pi() : 1.;

      // kinematic variables for PFJet associated to PFTau
      tauB->jetPt  = it->pfJetRef()->pt();
      tauB->jetEta = it->pfJetRef()->eta();
      tauB->jetPhi = it->pfJetRef()->phi();
    }
  }
  else {
    edm::LogError("TauBlock") << "Error! Can't get the product " << _inputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
