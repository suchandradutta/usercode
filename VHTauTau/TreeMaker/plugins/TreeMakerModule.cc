#include <cassert>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "VHTauTau/TreeMaker/plugins/TreeMakerModule.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "TTree.h"

TreeMakerModule::TreeMakerModule(const edm::ParameterSet& iConfig) : 
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _createTree(iConfig.getParameter<bool>("createTree"))
{
}
void TreeMakerModule::beginJob() 
{
  if (!_createTree) return;
  edm::Service<TFileService> fs;
  fs->file().cd("/");
  TTree* tree = fs->make<TTree>("vhtree", "");
  assert(tree);
  fs->file().ls();
}
//
// -- Analyze
//
void TreeMakerModule::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Get TTree pointer
  if (_createTree) return;
  TTree* tree = Utility::getTree("vhtree");
  tree->Fill();
}
void TreeMakerModule::endJob() {
  if (_createTree) return;
  TTree* tree = Utility::getTree("vhtree");
  tree->Write();  
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TreeMakerModule);
