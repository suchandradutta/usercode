#include "TTree.h"
#include "TClonesArray.h"
#include "VHTauTau/TreeMaker/plugins/GenJetBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

GenJetBlock::GenJetBlock(const edm::ParameterSet& iConfig) :
  _tree(0),
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("genJetSrc"))
{}
void GenJetBlock::beginJob() {
  // Get TTree pointer
  //edm::Service<TFileService> fs;
  //TTree* tree = fs->getObject<TTree>("vhtree", "treeCreator/");
  if (!_tree) _tree = Utility::getTree("vhtree");
  cloneGenJet = new TClonesArray("GenJet");
  _tree->Branch("GenJet", &cloneGenJet, 32000, 2);
  _tree->Branch("nGenJet", &fnGenJet,  "fnGenJet/I");
}
void GenJetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneGenJet->Clear();
  fnGenJet = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByLabel(_inputTag, genJets);

    if (genJets.isValid()) {
      edm::LogInfo("GenJetBlock") << "Total # GenJets: " << genJets->size();

      for (reco::GenJetCollection::const_iterator it = genJets->begin(); 
                                                 it != genJets->end(); ++it) {

        if (fnGenJet == kMaxGenJet) {
	  edm::LogInfo("GenJetBlock") << "Too many Gen Jets, fnGenJet = " << fnGenJet; 
	  break;
        }
        genJetB = new ((*cloneGenJet)[fnGenJet++]) GenJet();

        // fill in all the vectors
        genJetB->eta    = it->eta();
        genJetB->phi    = it->phi();
        genJetB->p      = it->p();
        genJetB->pt     = it->pt();
        genJetB->energy = it->energy();
        genJetB->emf    = it->emEnergy()/it->energy();
        genJetB->hadf   = it->hadEnergy()/it->energy();
      }
    } 
    else {
      edm::LogError("GenJetBlock") << "Error! Can't get the product " << _inputTag;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenJetBlock);
