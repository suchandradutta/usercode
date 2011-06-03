#include "TTree.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/METBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/MET.h"

METBlock::METBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("metSrc"))
{}
void METBlock::beginJob() 
{
  // Get TTree pointer
  //edm::Service<TFileService> fs;
  //TTree* tree = fs->getObject<TTree>("vhtree");
  TTree* tree = Utility::getTree("vhtree");
  cloneMET = new TClonesArray("MET");
  tree->Branch("MET", &cloneMET, 32000, 2);
  tree->Branch("nMET", &fnMET,  "fnMET/I");
}
void METBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneMET->Clear();
  fnMET = 0;

  edm::Handle<std::vector<pat::MET> > mets;
  iEvent.getByLabel(_inputTag, mets);

  if (mets.isValid()) {
    edm::LogInfo("METBlock") << "Total # PAT METs: " << mets->size();

    for (std::vector<pat::MET>::const_iterator it = mets->begin(); 
                                              it != mets->end(); ++it) {
      if (fnMET == kMaxMET) {
	edm::LogInfo("METBlock") << "Too many PAT MET, fnMET = " << fnMET; 
	break;
      }
      metB = new ((*cloneMET)[fnMET++]) MET();

      // fill in all the vectors
      metB->met          = it->pt();
      metB->metphi       = it->phi();
      metB->sumet        = it->sumEt();
      metB->metuncorr    = it->uncorrectedPt(pat::MET::uncorrALL);
      metB->metphiuncorr = it->uncorrectedPhi(pat::MET::uncorrALL);
      metB->sumetuncorr  = it->sumEt() - it->corSumEt(pat::MET::uncorrALL);
    }
  } 
  else {
    edm::LogError("METBlock") << "Error! Can't get the product " << _inputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METBlock);
