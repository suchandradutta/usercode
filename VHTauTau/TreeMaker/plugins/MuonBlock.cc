#include "TTree.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/MuonBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Math/GenVector/VectorUtil.h"

MuonBlock::MuonBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _muonInputTag(iConfig.getParameter<edm::InputTag>("muonSrc")),
  _pfMuonInputTag(iConfig.getParameter<edm::InputTag>("pfMuonSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  _beamSpotInputTag(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
  _beamSpotCorr(iConfig.getParameter<bool>("beamSpotCorr")),
  _muonID(iConfig.getParameter<std::string>("muonID"))
{}
void MuonBlock::beginJob() 
{
  // Get TTree pointer
  TTree* tree = Utility::getTree("vhtree");
  cloneMuon = new TClonesArray("Muon");
  tree->Branch("Muon", &cloneMuon, 32000, 2);
  tree->Branch("nMuon", &fnMuon,  "fnMuon/I");
}
void MuonBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneMuon->Clear();
  fnMuon = 0;

  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(_muonInputTag, muons);

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  //iEvent.getByLabel(_beamSpotInputTag, beamSpot);
  iEvent.getByLabel("offlineBeamSpot", beamSpot);

  edm::Handle<reco::PFCandidateCollection> pfMuons;
  iEvent.getByLabel(_pfMuonInputTag, pfMuons);

  //Define a Trigger Match Helper to retrieve trigger matches
  //  const pat::helper::TriggerMatchHelper matchHelper;
  
  if (muons.isValid()) {
    edm::LogInfo("MuonBlock") << "Total # Muons: " << muons->size();
    for (std::vector<pat::Muon>::const_iterator it = muons->begin(); it != muons->end(); ++it) {
      // if muon is not global muon, continue
      if (!it->isGlobalMuon()) continue;

      if (fnMuon == kMaxMuon) {
	edm::LogInfo("MuonBlock") << "Too many PAT Muons, fnMuon = " << fnMuon; 
	break;
      }
      double trkd0 = it->track()->d0();
      if (_beamSpotCorr && beamSpot.isValid()) trkd0 = -(it->track()->dxy(beamSpot->position()));
      else if (_beamSpotCorr && !beamSpot.isValid()) 
        edm::LogError("MuonsBlock") << "Error! Can't get the offlineBeamSpot";
      double reliso = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();

      // Vertex association
      double minVtxDist3D = 9999.;
      int indexVtx = -1;
      double vertexDistZ = 9999.;
      if (primaryVertices.isValid()) {
	edm::LogInfo("MuonsBlock") << "Total # Primary Vertices: " << primaryVertices->size();

        for (reco::VertexCollection::const_iterator v_it = primaryVertices->begin(); 
                                                   v_it != primaryVertices->end(); ++v_it) {
          double dist3D = sqrt(pow(it->track()->dxy(v_it->position()),2) 
                             + pow(it->track()->dz(v_it->position()),2));
          if (dist3D<minVtxDist3D) {
            minVtxDist3D = dist3D;
                indexVtx = int(std::distance(primaryVertices->begin(),v_it));
             vertexDistZ = it->track()->dz(v_it->position());
          }
        }
      } 
      else {
	edm::LogError("MuonBlock") << "Error! Can't get the product " << _vtxInputTag;
      }

      // PF-Isolation
      double pfreliso = it->userFloat("pfLooseIsoPt04")/it->pt();

      reco::TrackRef tk = it->innerTrack();

      muonB = new ((*cloneMuon)[fnMuon++]) Muon();
      muonB->eta        = it->eta();
      muonB->phi        = it->phi();
      muonB->pt         = it->pt();
      muonB->p          = it->p();
      muonB->energy     = it->energy();
      muonB->charge     = it->charge();
      muonB->trkD0      = trkd0;
      muonB->trkD0Error = it->track()->d0Error();
      muonB->trkDz      = it->track()->dz();
      muonB->trkDzError = it->track()->dzError();
      muonB->globalChi2 = it->track()->normalizedChi2();
      muonB->trkIso     = it->trackIso();
      muonB->ecalIso    = it->ecalIso();
      muonB->hcalIso    = it->hcalIso();
      muonB->hoIso      = it->isolationR03().hoEt;
      muonB->relIso     = reliso;
      muonB->passID     = (it->muonID(_muonID)) ? true : false;
      muonB->vtxDist3D  = minVtxDist3D;
      muonB->vtxIndex   = indexVtx;
      muonB->vtxDistZ   = vertexDistZ;
      const reco::HitPattern& hitp = tk->hitPattern(); 
      muonB->pixHits    = hitp.numberOfValidPixelHits();
      muonB->trkHits    = hitp.numberOfValidTrackerHits();
      muonB->muoHits    = hitp.numberOfValidMuonHits();
      muonB->matches    = it->numberOfMatches();
      muonB->pfRelIso   = pfreliso;

      muonB->isTrackerMuon = (it->isTrackerMuon()) ? true : false;

      // IP information
      muonB->dB  = it->dB(pat::Muon::PV2D);
      muonB->edB = it->edB(pat::Muon::PV2D);

      // UW recommendation
      muonB->isGlobalMuonPromptTight = muon::isGoodMuon(*it, muon::GlobalMuonPromptTight);
      muonB->isAllArbitrated         = muon::isGoodMuon(*it, muon::AllArbitrated);
      muonB->nChambers              = it->numberOfChambers();
      muonB->nMatches               = it->numberOfMatches();
      muonB->nMatchedStations       = it->numberOfMatchedStations();
      muonB->stationMask            = it->stationMask();
      muonB->stationGapMaskDistance = it->stationGapMaskDistance();
      muonB->stationGapMaskPull     = it->stationGapMaskPull();
    }
  } 
  else {
    edm::LogError("MuonBlock") << "Error! Can't get the product " << _muonInputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonBlock);
