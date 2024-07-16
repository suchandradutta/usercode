/********************************/
/*********************************/
/**                             **/
/**   Track Trigger Analysis    **/
/**                             **/
/**   suchandra.dutta@cern.ch   **/
/**                             **/
/*********************************/
/*********************************/
// system include files
#include <memory>
#include <cmath>
#include <typeinfo>
#include <vector>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>

// framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// user include files
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "L1Trigger/TrackTrigger/interface/TTClusterAlgorithm.h"
#include "L1Trigger/TrackTrigger/interface/TTClusterAlgorithmRecord.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Vertex.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"
#include "DataFormats/L1TCorrelator/interface/TkEm.h"
#include "DataFormats/L1TCorrelator/interface/TkEmFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkElectron.h"
#include "DataFormats/L1TCorrelator/interface/TkElectronFwd.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"

#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/TTAnalysis/interface/AnalObjects.h"

#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TMath.h"

#include <sstream>

using L1TTTrackType = TTTrack<Ref_Phase2TrackerDigi_>;
using L1TTTrackCollectionType = std::vector<L1TTTrackType>;
using L1TTStubCollection
      = std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>, TTStub<Ref_Phase2TrackerDigi_>>>;

// -----------------
// class declaration
// -----------------
class TTAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TTAnalysis(const edm::ParameterSet&);
  ~TTAnalysis();

private:
  void beginJob() override;
  //void beginRun(edm::Run const&, edm::EventSetup const&);
  void analyze(edm::Event const&, edm::EventSetup const&) override;
  void endJob() override = default;

  void fillGenParticleInfo();
  void fillSimTrackInfo();
  void fillL1TrackInfo();
  void fillL1EmuTrackInfo();
  void fillOfflineTrackInfo();
  void fillL1EGammaInfo(unsigned int label);

  void fillL1TkElectronInfo(unsigned int label, edm::Handle<l1t::TkElectronCollection>& handle);
  void fillL1TkEmInfo(unsigned int label, edm::Handle<l1t::TkEmCollection> & handle);
  void fillL1PFInfo();

  void fillL1PuppiInfo();
  bool matchGenInfo(const reco::Candidate* a, const reco::Candidate* b);
  int getSimTrkIndex(const unsigned int id);

  void findStubs(const TrackerTopology* topology, const L1TTTrackCollectionType::iterator& trk, TTStudy::Track& tkObj);
  // ----------member data ---------------------------
  // variables that fill ntuples:

  int nRecoTrack;
  TTStudy::Event* eventBr_;
  std::vector<TTStudy::Track>* l1TracksBr_;
  std::vector<TTStudy::Track>* l1EmuTracksBr_;
  std::vector<TTStudy::OfflineTrack>* offlineTracksBr_;
  std::vector<TTStudy::L1Object>* l1EGammaBr_;
  std::vector<TTStudy::L1TkObject>* l1tkEmuElecBr_;
  std::vector<TTStudy::L1TkObject>* l1tkEmuEmBr_;
  std::vector<TTStudy::L1PFObject>* l1pfCandBr_;
  std::vector<TTStudy::L1PFObject>* l1pupCandBr_;
  std::vector<TTStudy::GenParticle>* genParBr_;
  std::vector<TTStudy::SimTrack>* simTracksBr_;

  // ntuples:
  TTree* tree_;

  const edm::EDGetToken l1trackToken_;
  const edm::EDGetToken l1EmuTrackToken_;
  const edm::EDGetToken l1VertexToken_;
  const edm::EDGetToken l1EmuVertexToken_;
  std::vector<edm::EDGetToken> l1EGammaToken_;
  std::vector<edm::EDGetTokenT<l1t::TkElectronCollection>> l1tkEmuElecToken_;
  std::vector<edm::EDGetTokenT<l1t::TkEmCollection>> l1tkEmuEmToken_;
  const edm::EDGetToken l1pfCandToken_;
  const edm::EDGetToken l1pupCandToken_;
  const edm::EDGetToken offlineTrackToken_;
  const edm::EDGetToken offlineVertexToken_;
  const edm::EDGetToken beamSpotToken_;
  const edm::EDGetToken tpToken_;
  const edm::EDGetToken puToken_;
  const edm::EDGetToken simTrackToken_;
  const edm::EDGetToken simVertexToken_;
  const edm::EDGetToken genParticleToken_;

  const int debugFlag_;
  const bool l1TrackFlag_;
  const bool l1EmuTrackFlag_;
  const bool l1VertexFlag_;
  const bool l1EmuVertexFlag_;
  const bool l1EGammaFlag_;
  const bool l1tkEmuElecFlag_;
  const bool l1tkEmuEmFlag_;
  const bool l1pfCandFlag_;
  const bool l1pupCandFlag_;

  const bool offlineTrackFlag_;
  const bool offlineVertexFlag_;

  const bool simTrackFlag_;
  const bool genParticleFlag_;

  edm::ESHandle<MagneticField> theMagField;
  edm::ESHandle<TrackerGeometry> theTrkGeomHandle;

  edm::Handle<std::vector<l1t::Vertex>> l1VertexHandle_;
  edm::Handle<std::vector<l1t::VertexWord>> l1EmuVertexHandle_;
  edm::Handle<l1t::EGammaBxCollection> l1EGammaHandle_;
  edm::Handle<l1t::TkElectronCollection> l1tkEmuElecHandle_;
  edm::Handle<l1t::TkEmCollection> l1tkEmuEmHandle_;
  edm::Handle<l1t::PFCandidateCollection> l1pfCandHandle_;
  edm::Handle<l1t::PFCandidateCollection> l1pupCandHandle_;
  edm::Handle<L1TTTrackCollectionType> l1TTTrackHandle_;
  edm::Handle<L1TTTrackCollectionType> l1TTEmuTrackHandle_;
  edm::Handle<TrackingParticleCollection> tpHandle_;
  edm::Handle<reco::TrackCollection> offlineTrackHandle_;
  edm::Handle<reco::VertexCollection> offlineVertexHandle_;
  edm::Handle<reco::BeamSpot> beamSpotHandle_;
  edm::Handle<edm::SimTrackContainer> simTrackHandle_;
  edm::Handle<edm::SimVertexContainer> simVertexHandle_;
  edm::Handle<reco::GenParticleCollection> genParticleHandle_;

  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const TrackerTopology* tTopo_;

  reco::GenParticleCollection genParticles_;
  reco::TrackCollection offlineTracks_;
  edm::SimTrackContainer simTracks_;
  edm::SimVertexContainer simVertices_;
};

using std::cout;
using std::endl;
using std::cerr;
using std::vector;

// ----------------------------
// constructors and destructor
// ----------------------------
TTAnalysis::TTAnalysis(const edm::ParameterSet& iConfig) :
  l1trackToken_(consumes<vector<TTTrack<Ref_Phase2TrackerDigi_>>>(iConfig.getParameter<edm::InputTag>("L1TrackInputTag"))),
  l1EmuTrackToken_(consumes<vector<TTTrack<Ref_Phase2TrackerDigi_>>>(iConfig.getParameter<edm::InputTag>("L1EmulatedTrackInputTag"))),
  l1VertexToken_(consumes<vector<l1t::Vertex>>(iConfig.getParameter<edm::InputTag>("L1VertexInputTag"))),
  l1EmuVertexToken_(consumes<vector<l1t::VertexWord>>(iConfig.getParameter<edm::InputTag>("L1EmuVertexInputTag"))),
  l1pfCandToken_(consumes<vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("L1PFCandidateInputTag"))),
  l1pupCandToken_(consumes<vector<l1t::PFCandidate>>(iConfig.getParameter<edm::InputTag>("L1PuppiCandidateInputTag"))),
  offlineTrackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("OfflineTrackInputTag"))),
  offlineVertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("OfflineVertexInputTag"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotInputTag"))),
  tpToken_(consumes<vector<TrackingParticle>>(edm::InputTag("mix", "MergedTrackTruth"))),
  puToken_(consumes<vector<PileupSummaryInfo>>(edm::InputTag("addPileupInfo"))),
  simTrackToken_(consumes<edm::SimTrackContainer>(edm::InputTag("g4SimHits", ""))),
  simVertexToken_(consumes<edm::SimVertexContainer>(edm::InputTag("g4SimHits",""))),
  genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleInputTag"))),

  // flags
  debugFlag_       (iConfig.getParameter<int>("DebugFlag")),
  l1TrackFlag_     (iConfig.getParameter<bool>("L1TrackFlag")),
  l1EmuTrackFlag_  (iConfig.getParameter<bool>("L1EmulatedTrackFlag")),
  l1VertexFlag_    (iConfig.getParameter<bool>("L1VertexFlag")),
  l1EmuVertexFlag_ (iConfig.getParameter<bool>("L1EmuVertexFlag")),
  l1EGammaFlag_    (iConfig.getParameter<bool>("L1EGammaFlag")),
  l1tkEmuElecFlag_ (iConfig.getParameter<bool>("L1TkEmuElectronFlag")),
  l1tkEmuEmFlag_   (iConfig.getParameter<bool>("L1TkEmuEmFlag")),
  l1pfCandFlag_    (iConfig.getParameter<bool>("L1PFCandidateFlag")),
  l1pupCandFlag_    (iConfig.getParameter<bool>("L1PuppiCandidateFlag")),
  offlineTrackFlag_(iConfig.getParameter<bool>("OfflineTrackFlag")),
  offlineVertexFlag_(iConfig.getParameter<bool>("OfflineVertexFlag")),
  simTrackFlag_    (iConfig.getParameter<bool>("SimTrackFlag")),
  genParticleFlag_ (iConfig.getParameter<bool>("GenParticleFlag"))
{
  // now do whatever initialization is needed
  vector<edm::InputTag> egammaLabels = iConfig.getParameter<vector<edm::InputTag>>("L1EGammaInputTag");
  for (auto iTag = egammaLabels.begin(); iTag != egammaLabels.end(); ++iTag) {
    if (debugFlag_) cout << "--> Filling Input Tag for L1EGamma: " << *iTag << endl;
    l1EGammaToken_.emplace_back(consumes<l1t::EGammaBxCollection>(*iTag));
  }
  vector<edm::InputTag> l1tkEmuElecLabels = iConfig.getParameter<vector<edm::InputTag>>("L1TkEmuElectronInputTag");
  for (auto iTag = l1tkEmuElecLabels.begin(); iTag != l1tkEmuElecLabels.end(); ++iTag) {
    if (debugFlag_) cout << "--> Filling Input Tag for L1TkEmuElectron: " << *iTag << endl;
    l1tkEmuElecToken_.emplace_back(consumes<l1t::TkElectronCollection>(*iTag));
  }
  vector<edm::InputTag> l1tkEmuEmLabels = iConfig.getParameter<vector<edm::InputTag>>("L1TkEmuEmInputTag");
  for (auto iTag = l1tkEmuEmLabels.begin(); iTag != l1tkEmuEmLabels.end(); ++iTag) {
    if (debugFlag_) cout << "--> Filling Input Tag for L1TkEmuEm: " << *iTag << endl;
    l1tkEmuEmToken_.emplace_back(consumes<l1t::TkEmCollection>(*iTag));
  }

  topoToken_ = esConsumes<TrackerTopology, TrackerTopologyRcd>();

  usesResource(TFileService::kSharedResource);

  if (debugFlag_) {
    cout << "=> Analysis Flags: " << endl;
    cout << "\tReco Track Flag:             " << boolalpha << offlineTrackFlag_ << endl
	 << "\tReco Vertex Flag:            " << boolalpha << offlineVertexFlag_ << endl
         << "\tSim Track Flag:              " << boolalpha << simTrackFlag_ << endl
         << "\tGen Particle Flag:           " << boolalpha << genParticleFlag_ << endl
         << "\tL1 Track Flag:               " << boolalpha << l1TrackFlag_ << endl
         << "\tL1 Emulated Track Flag:      " << boolalpha << l1EmuTrackFlag_ << endl
         << "\tL1 Vertex Flag:              " << boolalpha << l1VertexFlag_ << endl
         << "\tL1 Tk Emulated Vertex Flag:  " << boolalpha << l1EmuVertexFlag_ << endl
         << "\tL1 EGamma Flag:              " << boolalpha << l1EGammaFlag_ << endl
         << "\tL1 Emulated TkElectron Flag: " << boolalpha << l1tkEmuElecFlag_ << endl
         << "\tL1 Emulated TkEm Flag:       " << boolalpha << l1tkEmuEmFlag_ << endl
	 << "\tL1 PFCadidate Flag:          " << boolalpha << l1pfCandFlag_ << endl
	 << "\tL1 PuppiCandidate Flag:      " << boolalpha << l1pupCandFlag_ << endl
	 << endl;
  }
}
TTAnalysis::~TTAnalysis() {
  delete eventBr_;
  if (l1TrackFlag_ && l1TracksBr_)           delete l1TracksBr_;
  if (l1EmuTrackFlag_ && l1EmuTracksBr_)     delete l1EmuTracksBr_;
  if (l1EGammaFlag_ && l1EGammaBr_)          delete l1EGammaBr_;
  if (l1tkEmuElecFlag_ && l1tkEmuElecBr_)    delete l1tkEmuElecBr_;
  if (l1tkEmuEmFlag_ && l1tkEmuEmBr_)        delete l1tkEmuEmBr_;
  if (l1pfCandFlag_ && l1pfCandBr_)          delete l1pfCandBr_;
  if (l1pupCandFlag_ && l1pupCandBr_)        delete l1pupCandBr_;
  if (offlineTrackFlag_ && offlineTracksBr_) delete offlineTracksBr_;
  if (genParticleFlag_ && genParBr_)         delete genParBr_;
  if (simTrackFlag_ && simTracksBr_)         delete simTracksBr_;
}
void TTAnalysis::beginJob() {
  //---------------------------------------------------------------------------------
  // Framework handles for the EVENTSETUP tracker geometry, L1 stack geometry, etc...
  // --------------------------------------------------------------------------------
  edm::Service<TFileService> fs;
  if (!fs.isAvailable()) return;

  tree_ = fs->make<TTree>("Phase2TriggerInfo", "Phase2 Trigger Information");

  eventBr_ = new TTStudy::Event();
  tree_->Branch("Event", "TTStudy::Event", &eventBr_); // , 32000, 2);

  if (l1TrackFlag_) {
    l1TracksBr_ = new vector<TTStudy::Track>();
    tree_->Branch("L1Track", "vector<TTStudy::Track>", &l1TracksBr_);
  }
  if (l1EmuTrackFlag_) {
    l1EmuTracksBr_ = new vector<TTStudy::Track>();
    tree_->Branch("L1EmulatedTrack", "vector<TTStudy::Track>", &l1EmuTracksBr_);
  }
  if (genParticleFlag_) {
    genParBr_ = new vector<TTStudy::GenParticle>();
    tree_->Branch("GenParticle", "vector<TTStudy::GenParticle>", &genParBr_);
  }
  if (simTrackFlag_) {
    simTracksBr_ = new vector<TTStudy::SimTrack>();
    tree_->Branch("SimTrack", "vector<TTStudy::SimTrack>", &simTracksBr_);
  }
  if (offlineTrackFlag_) {
    offlineTracksBr_ = new vector<TTStudy::OfflineTrack>();
    tree_->Branch("OfflineTrack", "vector<TTStudy::OfflineTrack>", &offlineTracksBr_);
  }
  if (l1EGammaFlag_) {
    l1EGammaBr_ = new vector<TTStudy::L1Object>();
    tree_->Branch("L1EGamma", "vector<TTStudy::L1Object>", &l1EGammaBr_);
  }
  if (l1tkEmuElecFlag_) {
    l1tkEmuElecBr_ = new vector<TTStudy::L1TkObject>();
    tree_->Branch("L1TkEmuElectron", "vector<TTStudy::L1TkObject>", &l1tkEmuElecBr_);
  }
  if (l1tkEmuEmFlag_) {
    l1tkEmuEmBr_ = new vector<TTStudy::L1TkObject>();
    tree_->Branch("L1TkEmuEm", "vector<TTStudy::L1TkObject>", &l1tkEmuEmBr_);
  }
  if (l1pfCandFlag_) {
    l1pfCandBr_ = new vector<TTStudy::L1PFObject>();
    tree_->Branch("L1PFCandidate", "std::vector<TTStudy::L1PFObject>", &l1pfCandBr_);
  }
  if (l1pupCandFlag_) {
    l1pupCandBr_ = new vector<TTStudy::L1PFObject>();
    tree_->Branch("L1PuppiCandidate", "std::vector<TTStudy::L1PFObject>", &l1pupCandBr_);
  }
}
void TTAnalysis::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  using namespace edm;

  // Tracker Topology
  tTopo_ = &iSetup.getData(topoToken_);

  // Clear the branches and get ready for the next event
  if (l1TrackFlag_ && l1TracksBr_)          l1TracksBr_->clear();
  if (offlineTrackFlag_&& offlineTracksBr_) offlineTracksBr_->clear();
  if (l1EmuTrackFlag_&& l1EmuTracksBr_)     l1EmuTracksBr_->clear();
  if (l1EGammaFlag_ && l1EGammaBr_)         l1EGammaBr_->clear();
  if (l1tkEmuElecFlag_ && l1tkEmuElecBr_)   l1tkEmuElecBr_->clear();
  if (l1tkEmuEmFlag_ && l1tkEmuEmBr_)       l1tkEmuEmBr_->clear();
  if (l1pfCandFlag_ && l1pfCandBr_)         l1pfCandBr_->clear();
  if (l1pupCandFlag_ && l1pupCandBr_)       l1pupCandBr_->clear();
  if (genParticleFlag_ && genParBr_)        genParBr_->clear();
  if (simTrackFlag_ && simTracksBr_)        simTracksBr_->clear();

  eventBr_->event = iEvent.id().event();
  //-----------------
  // Pile Up Vertices
  //-----------------
  int npu = 0;
  if (!iEvent.isRealData()) {
    edm::Handle<std::vector<PileupSummaryInfo>> puHandle;
    iEvent.getByToken(puToken_, puHandle);
    for (auto PVI = puHandle->begin(); PVI != puHandle->end(); ++PVI) {
      // More info about PU is here:
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#Accessing_PileupSummaryInfo_in_r
      if (debugFlag_ >> 4 & 0x1) cout << "bx = "  << PVI->getBunchCrossing() << ", nPU = " << PVI->getPU_NumInteractions() << endl;
      if (PVI->getBunchCrossing() == 0) npu = PVI->getPU_NumInteractions(); // in-time PU
      npu++;
    }
  }
  eventBr_->nPileUp = npu;

  //------------------------------
  // GenParticle/Truth information
  //------------------------------
  if (genParticleFlag_) {
    iEvent.getByToken(genParticleToken_, genParticleHandle_);
    if (genParticleHandle_.isValid()) {
      genParticles_ = (*genParticleHandle_.product());
      fillGenParticleInfo();
    }
    else {
      cerr << "analyze: GenParticleCollection for InputTag genParticles not found!" << endl;
    }
  }
  if (simTrackFlag_) {
    //--------------------------------
    // Get Sim Tracks and Sim Vertices
    //--------------------------------
    iEvent.getByToken(simTrackToken_, simTrackHandle_);
    if (simTrackHandle_.isValid()) {
      simTracks_ = (*simTrackHandle_.product());
      iEvent.getByToken(simVertexToken_, simVertexHandle_);
      if (simVertexHandle_.isValid()) {
	fillSimTrackInfo();
      }
      else {
	cerr << "analyze: SimVertex Collection for InputTag g4SimHits not found!" << endl;
      }
    }
    else {
      cerr << "analyze: SimTrack Collection for InputTag g4SimHits not found!" << endl;
    }
  }
  //-------------------
  // TrackingParticles
  //-------------------
  //  iEvent.getByLabel("mix", "MergedTrackTruth", tpHandle_);

  //------------------
  // Simulated L1Tracks
  //------------------
  if (l1TrackFlag_) {
    iEvent.getByToken(l1trackToken_, l1TTTrackHandle_);
    if (l1TTTrackHandle_.isValid()) fillL1TrackInfo();
    else {
      cerr << "analyze: L1 Track Collection not found!" << endl;
    }
  }
  //------------------
  // Emulated L1Tracks
  //------------------
  if (l1EmuTrackFlag_) {
    iEvent.getByToken(l1EmuTrackToken_, l1TTEmuTrackHandle_);
    if (l1TTEmuTrackHandle_.isValid()) fillL1EmuTrackInfo();
    else {
      cerr << "analyze: L1 Emulated Track Collection not found!" << endl;
    }
  }
  //-----------------------------------------
  // Offline Reconstructed Tracks and Vertex
  //-----------------------------------------
  if (offlineTrackFlag_) {
    iEvent.getByToken(offlineTrackToken_, offlineTrackHandle_);
    if (offlineTrackHandle_.isValid()) {
      iEvent.getByToken(offlineVertexToken_, offlineVertexHandle_);
      if (offlineVertexHandle_.isValid()) {
	const reco::Vertex& vit = offlineVertexHandle_->front();
	eventBr_->nOffVtx = (*offlineVertexHandle_.product()).size();
	eventBr_->zOffVtx = vit.z();
	eventBr_->sumPtOffVtx = vit.p4().P();
	fillOfflineTrackInfo();
      }
      else {
	cerr << "analyze: OfflinePrimaryVertex Collection not found!" << endl;
      }
    }
    else {
      cerr << "analyze: offline Track Collection not found!" << endl;
    }
  }
  //--------------------
  // L1 Primary Vertices
  //--------------------
  if (l1VertexFlag_) {
    edm::Handle<l1t::VertexCollection> handle;
    iEvent.getByToken(l1VertexToken_, handle);
    if (handle.isValid()) {
      vector<l1t::Vertex>::const_iterator it = handle->begin();
      eventBr_->zL1Vtx = it->z0();
      eventBr_->sumPtL1Vtx = it->pt();
      eventBr_->nL1Vtx = handle->size();
    }
    else {
      cerr << "analyze: l1t::Vertex Collection not found!" << endl;
    }
  }
  //---------------------
  // L1 Emulated Vertices
  //---------------------
  if (l1EmuVertexFlag_) {
    iEvent.getByToken(l1EmuVertexToken_, l1EmuVertexHandle_);
    if (l1EmuVertexHandle_.isValid()) {
      vector<l1t::VertexWord>::const_iterator it = l1EmuVertexHandle_->begin();
      eventBr_->zL1EmuVtx = it->z0Word();
      eventBr_->sumPtL1EmuVtx = it->ptWord();
      eventBr_->nL1EmuVtx     = l1EmuVertexHandle_->size();
    }
    else {
      cerr << "analyze: L1Tk Emulated Vertex Collection not found!" << endl;
    }
  }
  //----------
  // L1EGamma
  //----------
  if (l1EGammaFlag_) {
    for (auto i = 0U; i < l1EGammaToken_.size(); i++) {
      iEvent.getByToken(l1EGammaToken_[i], l1EGammaHandle_);
      if (l1EGammaHandle_.isValid()) fillL1EGammaInfo(i);
      else
	cerr << "analyze: L1EGamma Collection " << i << " not found!" << endl;
    }
  }
  //-----------------
  // L1TkEmuElectron
  //-----------------
  if (l1tkEmuElecFlag_) {
    for (auto i = 0U; i < l1tkEmuElecToken_.size(); i++) {
      iEvent.getByToken(l1tkEmuElecToken_[i], l1tkEmuElecHandle_);
      if (l1tkEmuElecHandle_.isValid()) fillL1TkElectronInfo(i,l1tkEmuElecHandle_);
      else
	cerr << "analyze: Emulated L1TkElectron Collection " << i << " not found!" << endl;
    }
  }
  //------------
  // L1TkEmuEm
  //------------
  if (l1tkEmuEmFlag_) {
    for (auto i = 0U; i < l1tkEmuEmToken_.size(); i++) {
      iEvent.getByToken(l1tkEmuEmToken_[i], l1tkEmuEmHandle_);
      if (l1tkEmuEmHandle_.isValid()) fillL1TkEmInfo(i,l1tkEmuEmHandle_);
      else
	cerr << "analyze: Emulated L1TkEm Collection " << i << " not found!" << endl;
    }
  }
  //----------------
  // L1PFCandidate
  //----------------
  if (l1pfCandFlag_) {
    iEvent.getByToken(l1pfCandToken_, l1pfCandHandle_);
    if (l1pfCandHandle_.isValid() ) fillL1PFInfo();
    else {
      cerr << "analyze: L1 PF Candidate not found!" << endl;
    }
  }
  //----------------
  // L1PuppiCandidate
  //----------------
  if (l1pupCandFlag_) {
    iEvent.getByToken(l1pupCandToken_, l1pupCandHandle_);
    if (l1pupCandHandle_.isValid() ) fillL1PuppiInfo();
    else {
      cerr << "analyze: L1 Puppi Candidate not found!" << endl;
    }
  }

  tree_->Fill();
}
int TTAnalysis::getSimTrkIndex(unsigned int id) {
  int index = -5;
  for (auto i = 0U; i < simTracks_.size(); i++) {
    auto simTrkId = simTracks_[i].trackId();
    if (simTrkId == id) {
      index = i;
      break;
    }
  }
  return index;
}
// ----------------------------------
// -- Fill GenParticle Information --
// ----------------------------------
void TTAnalysis::fillGenParticleInfo() {
  reco::GenParticleCollection::const_iterator gbeg = genParticles_.begin();
  if (debugFlag_ >> 2 & 0x1) {
    cout << setiosflags(ios::fixed)
         << setprecision(2)
         << "indx    status    pdgId  charge     eta      phi      pt     energy             mID                             dID"
	 << endl;
  }
  unsigned int i = 0;
  for (reco::GenParticleCollection::const_iterator it = genParticles_.begin(); it != genParticles_.end(); ++it) {
    TTStudy::GenParticle genPar;
    genPar.eta = it->eta();
    genPar.phi = it->phi();
    genPar.p = it->p();
    genPar.px = it->px();
    genPar.py = it->py();
    genPar.pz = it->pz();
    genPar.pt = it->pt();
    genPar.energy = it->energy();
    genPar.pdgId = it->pdgId();
    genPar.vx = it->vx();
    genPar.vy = it->vy();
    genPar.vz = it->vz();
    genPar.status = it->status();
    genPar.charge = it->charge();

    // First mother
    int idm = -1;
    const reco::Candidate* m = it->mother();
    if (m != nullptr) {
      for (reco::GenParticleCollection::const_iterator mit = genParticles_.begin(); mit != genParticles_.end(); ++mit) {
        const reco::Candidate* ap = &(*mit);
        if (matchGenInfo(m, ap)) {
          idm = distance(gbeg, mit);
          break;
        }
      }
    }
    genPar.motherIndex = idm;

    ostringstream dID;
    for (Size_t j = 0; j < it->numberOfDaughters(); ++j) {
      const reco::Candidate* d = it->daughter(j);
      for (reco::GenParticleCollection::const_iterator dit = genParticles_.begin(); dit != genParticles_.end(); ++dit) {
        const reco::Candidate* ap = &(*dit);
        if (matchGenInfo(d, ap)) {
          int idd = distance(gbeg, dit);
          genPar.daughterIndices.push_back(idd);
          dID << " " << idd;
          break;
        }
      }
    }
    if (debugFlag_ >> 2 & 0x1) {
      string ds = dID.str();
      if (!ds.length()) ds = " -";
      cout << setw(4)  << i++
	   << setw(8)  << it->status()
	   << setw(10) << it->pdgId()
	   << setw(8)  << it->charge()
	   << setw(10) << it->eta()
	   << setw(9)  << it->phi()
	   << setw(9)  << it->pt()
	   << setw(9)  << it->energy()
	   << setw(16) << idm
	   << ds
	   << endl;
    }
    genParBr_->push_back(genPar);
  }
  if (debugFlag_ >> 2 & 0x1) cout << resetiosflags(ios::fixed);
}
// ---------------------------------
// -- Fill SimTrack Information  ---
// ---------------------------------
void TTAnalysis::fillSimTrackInfo() {
  for (auto it = simTracks_.begin(); it != simTracks_.end(); ++it) {
    TTStudy::SimTrack simTk;
    simTk.pt = it->momentum().pt();
    simTk.eta = it->momentum().eta();
    simTk.phi = it->momentum().phi();
    int vertIndex = it->vertIndex();
    simTk.vtxIndx = vertIndex;
    if (static_cast<int>(simVertices_.size()) > vertIndex) {
      simTk.vx = simVertices_[vertIndex].position().x();
      simTk.vy = simVertices_[vertIndex].position().y();
      simTk.vz = simVertices_[vertIndex].position().z();
    }
    simTk.type = it->type();
    simTracksBr_->push_back(simTk);
  }
}
// --------------------------------------
// -- Fill Offline Track Information  ---
// --------------------------------------
void TTAnalysis::fillOfflineTrackInfo() {
  reco::TrackCollection offlineTracks = *(offlineTrackHandle_.product());
  for (auto iter = offlineTracks.begin(); iter != offlineTracks.end(); ++iter) {
    TTStudy::OfflineTrack trk;
    trk.pt  = iter->pt();
    trk.eta = iter->eta();
    trk.phi = iter->phi();
    trk.curvature = iter->qoverp();
    trk.chiSquare = iter->chi2(); // what about chi2Red
    trk.vertexX = iter->vx();
    trk.vertexY = iter->vy();
    trk.vertexZ = iter->vz();
    double trkd0 = iter->d0();
    double trkdz = iter->dz();
    if (beamSpotHandle_.isValid()) {
      trkd0 = -(iter->dxy(beamSpotHandle_->position()));
      trkdz = iter->dz(beamSpotHandle_->position());
    }
    trk.d0 = trkd0;
    trk.z0 = trkdz;

    double dxyWrtPV = 999;
    double dzWrtPV = 999;
    if (offlineVertexHandle_.isValid()) {
      const reco::Vertex& vit = offlineVertexHandle_->front();
      dxyWrtPV = iter->dxy(vit.position());
      dzWrtPV = iter->dz(vit.position());
    }
    trk.d0PV = dxyWrtPV;
    trk.z0PV = dzWrtPV;

    offlineTracksBr_->push_back(trk);
  }
}
bool TTAnalysis::matchGenInfo(const reco::Candidate* a, const reco::Candidate* b) {
  bool result = false;
  if ( a->pdgId()  == b->pdgId()  &&
       a->status() == b->status() &&
       a->pt()     == b->pt()     &&
       a->eta()    == b->eta()    &&
       a->phi()    == b->phi() ) result = true;
  return result;
}
// ---------------------------------
// -- Fill L1 Track Information  ---
// ---------------------------------
void TTAnalysis::fillL1TrackInfo() {
  L1TTTrackCollectionType l1TTTracks_ = *(l1TTTrackHandle_.product());
  if (debugFlag_ >> 1 & 0x1) cout << "Found " << l1TTTracks_.size() << " L1 Tracks" << endl;

  int itrk = 0;
  for (auto iter = l1TTTracks_.begin(); iter != l1TTTracks_.end() ; ++iter) {
    edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> l1track_ptr(l1TTTrackHandle_, itrk);
    itrk++;

    TTStudy::Track trk;
    trk.pt        = iter->momentum().perp();
    trk.eta       = iter->eta();
    trk.phi       = iter->phi();
    trk.localPhi  = iter->localPhi();
    trk.chiSquare = iter->chi2();
    trk.chiSquareRed = iter->chi2Red();
    trk.chi2BendRed  = iter->chi2BendRed();
    trk.curvature = iter->rInv();
    const GlobalPoint& gp = iter->POCA();
    trk.vertexX   = gp.x();
    trk.vertexY   = gp.y();
    trk.vertexZ   = gp.z();
    trk.d0        = iter->d0();
    trk.z0        = iter->z0();

    findStubs(tTopo_, iter, trk);

    trk.nFitPars = iter->nFitPars();
    trk.hitPattern = iter->hitPattern();
    l1TracksBr_->push_back(trk);
  }
}
// ------------------------------------------
// -- Fill L1 Emulated Track Information  ---
// ------------------------------------------
void TTAnalysis::fillL1EmuTrackInfo() {
  L1TTTrackCollectionType l1TTEmuTracks_ = *(l1TTEmuTrackHandle_.product());
  if (debugFlag_ >> 3 & 0x1) {
    cout << " Number of Tracks " << l1TTEmuTracks_.size() << endl
         << "  itrk       pT      Eta      Phi     chi2Bend  chi2RZ chi2RPhi"
         << " Curvature  VertexX  VertexY  VertexZ       d0       z0 nStub nStubPS"
	 << endl
         << setiosflags(ios::fixed)
         << setprecision(2);
  }
  int itrk = 0;
  for (auto iter = l1TTEmuTracks_.begin(); iter != l1TTEmuTracks_.end(); ++iter) {
    edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_>> l1track_ptr(l1TTEmuTrackHandle_, itrk);
    itrk++;
    TTStudy::Track trk;
    const GlobalVector& mom_trk = iter->momentum();  
    trk.pt        = mom_trk.perp();
    trk.eta       = mom_trk.eta();
    trk.phi       = mom_trk.phi();
    trk.localPhi  = iter->localPhi();
    trk.chiSquare = iter->chi2();
    trk.chiSquareRed = iter->chi2Red();
    trk.chi2BendRed  = iter->stubPtConsistency(); //iter->chi2BendRed();
    trk.chi2RZRed    = iter->chi2ZRed();
    trk.chi2RPhiRed  = iter->chi2XYRed();

    trk.curvature = iter->rInv();

    const GlobalPoint& gp = iter->POCA();
    trk.vertexX   = gp.x();
    trk.vertexY   = gp.y();
    trk.vertexZ   = gp.z();

    trk.d0        = iter->d0();
    trk.z0        = iter->z0();

    findStubs(tTopo_, iter, trk);

    trk.mva1 = iter->trkMVA1();
    trk.mva2 = iter->trkMVA2();
    trk.mva3 = iter->trkMVA3();

    trk.nFitPars = iter->nFitPars();
    trk.hitPattern = iter->hitPattern();

    l1EmuTracksBr_->push_back(trk);

    if (debugFlag_ >> 3 & 0x1) {
      //cout << "nPar: " << iter->nFitPars() << ", d0: " << iter->d0() << endl;
      cout << setw(6) << itrk
	   << setprecision(2)
	   << setw(9) << trk.pt
	   << setw(9) << trk.eta
	   << setw(9) << trk.phi
	   << setw(9) << trk.chi2BendRed
	   << setprecision(3)
	   << setw(9) << trk.chi2RZRed
	   << setw(9) << trk.chi2RPhiRed
	   << setprecision(4)
	   << setw(10) << trk.curvature
	   << setprecision(3)
	   << setw(9) << trk.vertexX
	   << setw(9) << trk.vertexY
	   << setw(9) << trk.vertexZ
	   << setw(9) << trk.d0
	   << setw(9) << trk.z0
	   << setw(6) << trk.nStub
	   << setw(8) << trk.nStub_PS
	   << endl;
    }
  }
  if (debugFlag_ >> 3 & 0x1) cout << resetiosflags(ios::fixed);
}
void TTAnalysis::findStubs(const TrackerTopology* topology,
			   const L1TTTrackCollectionType::iterator& trk,
			   TTStudy::Track& tkObj)
{
  const auto& stubs = trk->getStubRefs();
  tkObj.nStub = stubs.size();

  int nStub_PS = 0, nStub_SS = 0;
  for (auto it  = stubs.begin(); it != stubs.end(); ++it) {
    auto detid = (*it)->getDetId();
    if (detid.det() != DetId::Detector::Tracker) continue;
    if (detid.subdetId() == StripSubdetector::TOB) {
      (topology->tobLayer(detid) <= 3) ? nStub_PS++ : nStub_SS++;
    }
    else if (detid.subdetId() == StripSubdetector::TID) {
      (topology->tidRing(detid) <= 9) ? nStub_PS++ : nStub_SS++;
    }
  }
  tkObj.nStub_PS = nStub_PS;
  tkObj.nStub_SS = nStub_SS;
}
void TTAnalysis::fillL1EGammaInfo(unsigned int label) {
  if (debugFlag_ >> 1 & 0x1) cout << " -----  Accessing L1EGamma objects ---- of type " << label << endl;
  l1t::EGammaBxCollection eGammaCollection = (*l1EGammaHandle_.product());
  int i = 0;
  for (auto iter = eGammaCollection.begin(0); iter != eGammaCollection.end(0);  ++iter) { // considering BX = only
    i++;
    TTStudy::L1Object l1EGamma;
    l1EGamma.e      = iter->energy();
    l1EGamma.et     = iter->et();
    l1EGamma.phi    = iter->phi();
    l1EGamma.eta    = iter->eta();
    l1EGamma.label  = label;
    l1EGamma.hwQual = iter->hwQual();

    l1EGammaBr_->push_back(l1EGamma);
  }
}
void TTAnalysis::fillL1TkElectronInfo(unsigned int label,
				      edm::Handle<l1t::TkElectronCollection>& handle)
{
  if (debugFlag_ >> 1 & 0x1) cout << " -----  Accessing L1TkElectronPaticle objects ---- " << endl;
  for (auto iter = handle->begin(); iter != handle->end(); ++iter) {
    TTStudy::L1TkObject l1tkElec;
    l1tkElec.e = iter->energy();
    l1tkElec.et = iter->et();
    l1tkElec.phi = iter->phi();
    l1tkElec.eta = iter->eta();
    l1tkElec.z0 = iter->trkzVtx();

    const edm::Ptr<L1TTTrackType>& trk = iter->trkPtr();
    const GlobalVector& mom_trk = trk->momentum();  
    l1tkElec.pt_tk  = mom_trk.perp();
    l1tkElec.eta_tk = mom_trk.eta();
    l1tkElec.phi_tk = mom_trk.phi();
    l1tkElec.chiSquare = trk->chi2();
    l1tkElec.chiSquareRed = trk->chi2Red();
    l1tkElec.curvature = trk->rInv();

    const GlobalPoint& gp = trk->POCA();
    l1tkElec.vertexX   = gp.x();
    l1tkElec.vertexY   = gp.y();
    l1tkElec.d0        = gp.x() * sin(gp.phi()) + gp.y() * cos(gp.phi());

    l1tkElec.trackIsolation = iter->trkIsol();
    l1tkElec.trackIsolationPV = iter->trkIsolPV();
    l1tkElec.pfIsolation = iter->pfIsol();
    l1tkElec.pfIsolationPV = iter->pfIsolPV();
    l1tkElec.puppiIsolation = iter->puppiIsol();
    l1tkElec.puppiIsolationPV = iter->puppiIsolPV();
    l1tkElec.hwQual = iter->egCaloPtr()->hwQual();
    l1tkElec.label  = label;

    l1tkEmuElecBr_->push_back(l1tkElec);
  }
}
void TTAnalysis::fillL1TkEmInfo(unsigned int label, edm::Handle<l1t::TkEmCollection>& handle) {
  if (debugFlag_ >> 1 & 0x1) cout << " -----  Accessing L1TkEmPaticle objects ---- " << endl;
  for (auto iter = handle->begin(); iter != handle->end(); ++iter) {
    TTStudy::L1TkObject l1tkEm;
    l1tkEm.e = iter->energy();
    l1tkEm.et = iter->et();
    l1tkEm.phi = iter->phi();
    l1tkEm.eta = iter->eta();
    l1tkEm.z0 = -1.0;

    l1tkEm.pt_tk  = -99.9;
    l1tkEm.eta_tk = -99.9;
    l1tkEm.phi_tk = -99.9;
    l1tkEm.chiSquare = -99.9;
    l1tkEm.chiSquareRed = -99.9;
    l1tkEm.curvature = -99.9;
    l1tkEm.vertexX   = -99.9;
    l1tkEm.vertexY   = -99.9;
    l1tkEm.d0        = -99.9;

    l1tkEm.trackIsolation   = iter->trkIsol();
    l1tkEm.trackIsolationPV = iter->trkIsolPV();
    l1tkEm.pfIsolation      = iter->pfIsol();
    l1tkEm.pfIsolationPV    = iter->pfIsolPV();
    l1tkEm.puppiIsolation   = iter->puppiIsol();
    l1tkEm.puppiIsolationPV = iter->puppiIsolPV();

    l1tkEm.hwQual = iter->egCaloPtr()->hwQual();
    l1tkEm.label  = label;

    l1tkEmuEmBr_->push_back(l1tkEm);
  }
}
void TTAnalysis::fillL1PFInfo() {
  if (debugFlag_ >> 1 & 0x1) cout << " -----  Accessing L1PFCanidate objects ---- " << endl;
  const l1t::PFCandidateCollection& candCollection = *(l1pfCandHandle_.product());
  for (const auto& cand: candCollection) {
    TTStudy::L1PFObject l1Obj;
    l1Obj.e   = cand.energy();
    l1Obj.pt  = cand.pt();
    l1Obj.eta = cand.eta();
    l1Obj.phi = cand.phi();
    l1Obj.charge = cand.charge();
    l1Obj.dxy = cand.dxy();
    l1Obj.z0  = cand.z0();
    l1Obj.puppyWeight = cand.puppiWeight();
    l1Obj.pdgId = cand.pdgId();
    if (cand.pfTrack().isNonnull()) {
      l1Obj.pt_refTrk   = cand.pfTrack()->pt();
      l1Obj.eta_refTrk  = cand.pfTrack()->eta();
      l1Obj.phi_refTrk  = cand.pfTrack()->phi();
    }
    if (cand.pfCluster().isNonnull()) {
      l1Obj.pt_refClus  = cand.pfCluster()->pt();
      l1Obj.eta_refClus = cand.pfCluster()->eta();
      l1Obj.phi_refClus = cand.pfCluster()->phi();
    }
    l1pfCandBr_->push_back(l1Obj);
  }
}
void TTAnalysis::fillL1PuppiInfo() {
  if (debugFlag_ >> 1 & 0x1) cout << " -----  Accessing L1PFCanidate objects ---- " << endl;
  const l1t::PFCandidateCollection& candCollection = *(l1pupCandHandle_.product());
  for (const auto& cand:candCollection) {
    TTStudy::L1PFObject l1Obj;
    l1Obj.e   = cand.energy();
    l1Obj.pt  = cand.pt();
    l1Obj.eta = cand.eta();
    l1Obj.phi = cand.phi();
    l1Obj.charge = cand.charge();
    l1Obj.dxy = cand.dxy();
    l1Obj.z0  = cand.z0();
    l1Obj.puppyWeight = cand.puppiWeight();
    l1Obj.pdgId = cand.pdgId();
    if (cand.pfTrack().isNonnull()) {
      l1Obj.pt_refTrk   = cand.pfTrack()->pt();
      l1Obj.eta_refTrk  = cand.pfTrack()->eta();
      l1Obj.phi_refTrk  = cand.pfTrack()->phi();
    }
    if (cand.pfCluster().isNonnull()) {
      l1Obj.pt_refClus  = cand.pfCluster()->pt();
      l1Obj.eta_refClus = cand.pfCluster()->eta();
      l1Obj.phi_refClus = cand.pfCluster()->phi();
    }
    l1pupCandBr_->push_back(l1Obj);
  }
}
// define this as a plug-in
DEFINE_FWK_MODULE(TTAnalysis);
