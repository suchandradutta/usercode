#include "TTree.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/GenMETBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"

GenMETBlock::GenMETBlock(const edm::ParameterSet& iConfig) :
  _tree(0),
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("genMETSrc"))
{}
void GenMETBlock::beginJob() 
{
  // Get TTree pointer
  //edm::Service<TFileService> fs;
  //TTree* tree = fs->getObject<TTree>("vhtree");

  if (!_tree) _tree = Utility::getTree("vhtree");
  cloneGenMET = new TClonesArray("GenMET");
  _tree->Branch("GenMET", &cloneGenMET, 32000, 2);
  _tree->Branch("nGenMET", &fnGenMET,  "fnGenMET/I");
}
void GenMETBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneGenMET->Clear();
  fnGenMET = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenMETCollection> mets;
    iEvent.getByLabel(_inputTag, mets);
    if (mets.isValid()) {
      edm::LogInfo("GenMETBlock") << "Total # GenMETs: " << mets->size();
      for (reco::GenMETCollection::const_iterator it = mets->begin(); 
                                                 it != mets->end(); ++it ) {
        if (fnGenMET == kMaxGenMET) {
	  edm::LogInfo("GenMETBlock") << "Too many GenMET, fnGenMET = " << fnGenMET; 
	  break;
        }
        genMetB = new ((*cloneGenMET)[fnGenMET++]) GenMET();

        // fill in all the vectors
        genMetB->met    = it->pt();
        genMetB->metphi = it->phi();
        genMetB->sumet  = it->sumEt();
      }
    } 
    else {
      edm::LogError("GenMETBlock") << "Error! Can't get the product " << _inputTag;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenMETBlock);
