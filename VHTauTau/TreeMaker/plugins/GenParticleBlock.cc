#include "TTree.h"
#include "TClonesArray.h"
#include "VHTauTau/TreeMaker/plugins/GenParticleBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

GenParticleBlock::GenParticleBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("genParticleSrc"))
{}
void GenParticleBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = Utility::getTree("vhtree");
  cloneGenParticle = new TClonesArray("GenParticle");
  tree->Branch("GenParticle", &cloneGenParticle, 32000, 2);
  tree->Branch("nGenParticle", &fnGenParticle,  "fnGenParticle/I");
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
        genParticleB->status    = it->status();
        genParticleB->numDaught = it->numberOfDaughters();
        genParticleB->numMother = it->numberOfMothers();

        int idx = -1;
        for (reco::GenParticleCollection::const_iterator mit = genParticles->begin(); 
                                                        mit != genParticles->end(); ++mit ) {
          if ( it->mother() == &(*mit) ) {
	    idx = std::distance(genParticles->begin(),mit);
	    break;
          }
        }
        genParticleB->motherIndex = idx;

	if (_verbosity > 0) std::cout << "pdgId=" << it->pdgId() << std::endl;
        for (size_t j = 0; j < it->numberOfMothers(); ++j) {
	  const reco::Candidate* m = it->mother(j);
          for (reco::GenParticleCollection::const_iterator mit = genParticles->begin(); 
                                                          mit != genParticles->end(); ++mit) {
            if (m == &(*mit) ) { // -- keep all the entries && it->pdgId() != mit->pdgId() ) {
	      int idx = std::distance(genParticles->begin(), mit);
	      if (_verbosity > 0) std::cout << "mother index/pdgId: " << idx << "/" << mit->pdgId() << std::endl;
              genParticleB->motherIndices.push_back(idx);
              break;
            }
          } 
        }
        for (size_t j = 0; j < it->numberOfDaughters(); ++j) {
	  const reco::Candidate* d = it->daughter(j);
          for (reco::GenParticleCollection::const_iterator mit = genParticles->begin(); 
                                                          mit != genParticles->end(); ++mit) {
            if (d == &(*mit) ) { // -- keep all the entries  && it->pdgId() != mit->pdgId() ) {
	      int idx = std::distance(genParticles->begin(), mit);
	      if (_verbosity > 0) std::cout << "daughter index/pdgId: " << idx << "/" << mit->pdgId() << std::endl;
              genParticleB->daughtIndices.push_back(idx);
              break;
            }
          } 
        }
        if (_verbosity > 0) std::cout << "# of mother/daughter: "  
                  << genParticleB->motherIndices.size() << "/" 
		  << genParticleB->daughtIndices.size() << std::endl;
      }
    } 
    else {
      edm::LogError("GenParticleBlock") << "Error! Can't get the product " << _inputTag;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenParticleBlock);
