#include <iostream>

#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Math/GenVector/VectorUtil.h"

#include "VHTauTau/TreeMaker/plugins/MuonBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

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
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneMuon = new TClonesArray("vhtm::Muon");
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
  iEvent.getByLabel(_beamSpotInputTag, beamSpot);

  edm::Handle<reco::PFCandidateCollection> pfMuons;
  iEvent.getByLabel(_pfMuonInputTag, pfMuons);

  if (muons.isValid()) {
    edm::LogInfo("MuonBlock") << "Total # Muons: " << muons->size();
    for (std::vector<pat::Muon>::const_iterator it  = muons->begin(); 
                                                it != muons->end(); ++it) {
      // if it is not a global muon, continue
      if (!it->isGlobalMuon()) continue;

      if (fnMuon == kMaxMuon) {
	edm::LogInfo("MuonBlock") << "Too many PAT Muons, fnMuon = " << fnMuon; 
	break;
      }

      reco::TrackRef tk  = it->innerTrack();  // tracker segment only
      reco::TrackRef gtk = it->globalTrack(); 

      double trkd0 = gtk->d0();
      double trkdz = gtk->dz();
      if (_beamSpotCorr) {
        if (beamSpot.isValid()) {
          trkd0 = -(gtk->dxy(beamSpot->position()));
          trkdz = gtk->dz(beamSpot->position());
        }
        else
          edm::LogError("MuonsBlock") << "Error >> Failed to get BeamSpot for label: "
                                      << _beamSpotInputTag;
      }
      double reliso = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();

      // Vertex association
      double minVtxDist3D = 9999.;
      int indexVtx = -1;
      double vertexDistZ = 9999.;
      if (primaryVertices.isValid()) {
	edm::LogInfo("MuonBlock") << "Total # Primary Vertices: " << primaryVertices->size();

        for (reco::VertexCollection::const_iterator v_it  = primaryVertices->begin(); 
                                                    v_it != primaryVertices->end(); ++v_it) {
          double dxy = gtk->dxy(v_it->position());
          double dz  = gtk->dz(v_it->position());
          double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
          if (dist3D < minVtxDist3D) {
            minVtxDist3D = dist3D;
            indexVtx = int(std::distance(primaryVertices->begin(), v_it));
            vertexDistZ = dz;
          }
        }
      } 
      else {
	edm::LogError("MuonBlock") << "Error >> Failed to get VertexCollection for label: " 
                                   << _vtxInputTag;
      }

      // PF-Isolation
      double pfreliso = it->userFloat("pfLooseIsoPt04")/it->pt();

      // UW prescription for pf based isolation
      //(mu.chargedHadronIso" "+max(mu.photonIso()" "+mu.neutralHadronIso()" "-0.5*mu.userIso(0),0.0))" "/mu.pt() < ${threshold}"
      double v1 = it->photonIso() + it->neutralHadronIso() - 0.5*it->userIso(0);
      double vmax = (v1 > 0) ? v1 : 0;
      double UWpfreliso = (it->chargedHadronIso() + vmax)/it->pt();

      muonB = new ((*cloneMuon)[fnMuon++]) vhtm::Muon();
      muonB->eta        = it->eta();
      muonB->phi        = it->phi();
      muonB->pt         = it->pt();
      muonB->p          = it->p();
      //muonB->ptError    = gtk->ptError();
      muonB->energy     = it->energy();
      muonB->charge     = it->charge();
      muonB->trkD0      = trkd0;
      muonB->trkD0Error = gtk->d0Error();
      muonB->trkDz      = trkdz;
      muonB->trkDzError = gtk->dzError();
      muonB->globalChi2 = it->normChi2();
      muonB->trkIso     = it->trackIso();
      muonB->ecalIso    = it->ecalIso();
      muonB->hcalIso    = it->hcalIso();
      muonB->hoIso      = it->isolationR03().hoEt;
      muonB->relIso     = reliso;
      muonB->passID     = (it->muonID(_muonID)) ? true : false;
      muonB->vtxDist3D  = minVtxDist3D;
      muonB->vtxIndex   = indexVtx;
      muonB->vtxDistZ   = vertexDistZ;
      const reco::HitPattern& hitp = gtk->hitPattern(); 
      muonB->pixHits    = hitp.numberOfValidPixelHits();
      muonB->trkHits    = hitp.numberOfValidTrackerHits();
      muonB->muoHits    = hitp.numberOfValidMuonHits();
      muonB->matches    = it->numberOfMatches();
      muonB->pfRelIso   = UWpfreliso;

      muonB->isTrackerMuon = (it->isTrackerMuon()) ? true : false;

      // IP information
      muonB->dB  = it->dB(pat::Muon::PV2D);
      muonB->edB = it->edB(pat::Muon::PV2D);

      muonB->dB3d  = it->dB(pat::Muon::PV3D);
      muonB->edB3d = it->edB(pat::Muon::PV3D);

      // UW recommendation
      muonB->isGlobalMuonPromptTight = muon::isGoodMuon(*it, muon::GlobalMuonPromptTight);
      muonB->isAllArbitrated         = muon::isGoodMuon(*it, muon::AllArbitrated);
      muonB->nChambers              = it->numberOfChambers();
      muonB->nMatches               = it->numberOfMatches();
      muonB->nMatchedStations       = it->numberOfMatchedStations();
      muonB->stationMask            = it->stationMask();
      muonB->stationGapMaskDistance = it->stationGapMaskDistance();
      muonB->stationGapMaskPull     = it->stationGapMaskPull();

      // Vertex information
      //const reco::Candidate::Point& vertex = it->vertex();
      //muonB->vx = vertex.x();             
      //muonB->vy = vertex.y();             
      //muonB->vz = vertex.z();             
    }
  } 
  else {
    edm::LogError("MuonBlock") << "Error >> Failed to get pat::Muon collection for label: " 
                               << _muonInputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonBlock);
