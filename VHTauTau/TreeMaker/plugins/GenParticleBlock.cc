#include "TTree.h"
#include "TClonesArray.h"
#include "VHTauTau/TreeMaker/plugins/GenParticleBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

GenParticleBlock::GenParticleBlock(const edm::ParameterSet& iConfig) :
  _tree(0),
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("genParticleSrc"))
{}
void GenParticleBlock::beginJob() {
  // Get TTree pointer
  //  edm::Service<TFileService> fs;
  //TTree* tree = fs->getObject<TTree>("vhtree");
  if (!_tree) _tree = Utility::getTree("vhtree");
  cloneGenParticle = new TClonesArray("GenParticle");
  _tree->Branch("GenParticle", &cloneGenParticle, 32000, 2);
  _tree->Branch("nGenParticle", &fnGenParticle,  "fnGenParticle/I");
}
void GenParticleBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneGenParticle->Clear();
  fnGenParticle = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(_inputTag, genParticles);
    if (genParticles.isValid()) {
      edm::LogInfo("GenParticleBlock") << "Total # GenParticles: " << genParticles->size();
      for (reco::GenParticleCollection::const_iterator it = genParticles->begin(); 
                                                      it != genParticles->end(); ++it ) {
        if (fnGenParticle == kMaxGenParticle) {
	  edm::LogInfo("GenParticleBlock") << "Too many GenParticles, fnGenParticle = " << fnGenParticle;
	  break;
        }
        genParticleB = new ((*cloneGenParticle)[fnGenParticle++]) GenParticle();

        // fill in all the vectors
        genParticleB->eta       = it->eta();
        genParticleB->phi       = it->phi();
        genParticleB->p         = it->p();
        genParticleB->px        = it->px();
        genParticleB->py        = it->py();
        genParticleB->pz        = it->pz();
        genParticleB->pt        = it->pt();
        genParticleB->energy    = it->energy();
        genParticleB->pdgId     = it->pdgId();
        genParticleB->vx        = it->vx();
        genParticleB->vy        = it->vy();
        genParticleB->vz        = it->vz();
        genParticleB->numDaught = it->numberOfDaughters();
        genParticleB->status    = it->status();

        int idx = -1;
        for (reco::GenParticleCollection::const_iterator mit = genParticles->begin(); 
                                                        mit != genParticles->end(); ++mit) {
          if (it->mother() == &(*mit)) {
	    idx = std::distance(genParticles->begin(), mit);
	    break;
          }
        }
        genParticleB->motherIndex = idx;
      }
    } 
    else {
      edm::LogError("GenParticleBlock") << "Error! Can't get the product " << _inputTag;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenParticleBlock);
