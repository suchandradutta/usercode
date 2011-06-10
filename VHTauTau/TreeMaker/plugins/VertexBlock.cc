#include "TTree.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/VertexBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

VertexBlock::VertexBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("vertexSrc"))
{}
void VertexBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = Utility::getTree("vhtree");
  cloneVertex = new TClonesArray("Vertex");
  tree->Branch("Vertex", &cloneVertex, 32000, 2);
  tree->Branch("nVertex", &fnVertex, "fnVertex/I");
}
void VertexBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneVertex->Clear();
  fnVertex = 0;

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_inputTag, primaryVertices);

  if (primaryVertices.isValid()) {
    edm::LogInfo("VertexBlock") << "Total # Primary Vertices: " << primaryVertices->size();

    for (reco::VertexCollection::const_iterator it = primaryVertices->begin(); 
                                               it != primaryVertices->end(); ++it) {
      if (fnVertex == kMaxVertex) {
	edm::LogInfo("VertexBlock") << "Too many Vertex, fnVertex = " << fnVertex; 
	break;
      }
      vertexB = new ((*cloneVertex)[fnVertex++]) Vertex();
      vertexB->x = it->x();
      vertexB->y = it->y();
      vertexB->z = it->z();
      vertexB->xErr = it->xError();
      vertexB->yErr = it->yError();
      vertexB->zErr = it->zError();
      vertexB->rho =  it->position().rho();
      vertexB->chi2 = it->chi2();
      vertexB->ndf = it->ndof();
      vertexB->ntracks = int(it->tracksSize());
      vertexB->ntracksw05 = it->nTracks(0.5); // number of tracks in the vertex with weight above 0.5
      vertexB->isfake = it->isFake();
    }
  } 
  else {
    edm::LogError("VertexBlock") << "Error! Can't get the product " << _inputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VertexBlock);
