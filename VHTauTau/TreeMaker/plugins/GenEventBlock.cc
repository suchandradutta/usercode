#include <algorithm>
#include <iostream>

#include "TTree.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/GenEventBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

GenEventBlock::GenEventBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  _genEvtInfoInputTag(iConfig.getParameter<edm::InputTag>("GenEventInfoInputTag")),
  _storePDFWeights(iConfig.getParameter<bool>("StorePDFWeights")),
  _pdfWeightsInputTag(iConfig.getParameter<edm::InputTag>("PDFWeightsInputTag"))
{}
GenEventBlock::~GenEventBlock() {
  delete _pdfWeights;
}
void GenEventBlock::beginJob() 
{
  _pdfWeights = new std::vector<double>();

  // Get TTree pointer
  TTree* tree = Utility::getTree("vhtree");
  cloneGenEvent = new TClonesArray("GenEvent");
  tree->Branch("GenEvent", &cloneGenEvent, 32000, 2);
  tree->Branch("pdfWeights", "vector<double>", &_pdfWeights);
}
void GenEventBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray 
  cloneGenEvent->Clear();

  // Clear the two independent vectors
  _pdfWeights->clear();

  if (!iEvent.isRealData()) {
    // GenEventInfo Part
    edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
    iEvent.getByLabel(_genEvtInfoInputTag, genEvtInfoProduct);
    if (genEvtInfoProduct.isValid()) {
      edm::LogInfo("GenEventBlock") << "Successfully obtained " << _genEvtInfoInputTag;
      genEventB->processID = genEvtInfoProduct->signalProcessID();
      genEventB->ptHat     = (genEvtInfoProduct->hasBinningValues() ? genEvtInfoProduct->binningValues()[0] : 0.);
    } 
    else {
      edm::LogError("GenEventBlock") << "Error! Can't get the product " << _genEvtInfoInputTag;
    }
    // PDF Weights Part
    if (_storePDFWeights) {
      edm::Handle<std::vector<double> > pdfWeightsHandle;
      iEvent.getByLabel(_pdfWeightsInputTag, pdfWeightsHandle);
      if (pdfWeightsHandle.isValid()) {
	edm::LogInfo("GenEventBlock") << "Successfully obtained " << _pdfWeightsInputTag;
	copy(pdfWeightsHandle->begin(), pdfWeightsHandle->end(), genEventB->pdfWeights.begin());
	copy(pdfWeightsHandle->begin(), pdfWeightsHandle->end(), _pdfWeights->begin());
      } 
      else {
	edm::LogError("GenEventBlock") << "Error! Can't get the product " << _pdfWeightsInputTag;
      }
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenEventBlock);
