#include "TTree.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/EventSkimmer.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/Common/interface/Handle.h"

#include "Math/GenVector/VectorUtil.h"

EventSkimmer::EventSkimmer(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _muonInputTag(iConfig.getParameter<edm::InputTag>("muonSrc")),
  _electronInputTag(iConfig.getParameter<edm::InputTag>("electronSrc")),
  _doMuonSelection(iConfig.getParameter<bool>("doMuonSelection")),
  _doElectronSelection(iConfig.getParameter<bool>("doElectronSelection")),
  _muonPtCut(iConfig.getParameter<double>("muonPtCut")),
  _electronPtCut(iConfig.getParameter<double>("electronPtCut")),
  _createHisto(iConfig.getParameter<bool>("createHisto"))
{}
void EventSkimmer::beginJob() 
{
  _totalEvents    = 0;
  _selectedEvents = 0;
}
bool EventSkimmer::filter(edm::Event& iEvent, edm::EventSetup const& iSetup) {
  _totalEvents++;
  bool retVal = false; 
  if (_doMuonSelection)          retVal = checkMuonSelection(iEvent);
  else if (_doElectronSelection) retVal = checkElectronSelection(iEvent);
  if (retVal) _selectedEvents++; 
  return retVal;
}
void EventSkimmer::endJob() {
  if (_createHisto) {
    edm::Service<TFileService> fs;
    fs->file().cd("/");
    TH1F* th = fs->make<TH1F>("EventSkimInfo", "Information about Event Skim " , 2, -0.5, 1.5);
    th->SetBinContent(1, _totalEvents);
    th->SetBinContent(2, _selectedEvents);
  }
  edm::LogInfo("EventSkimmer") << " Total # of Events "<< _totalEvents;
  edm::LogInfo("EventSkimmer") << " Selected # Events "<< _selectedEvents;
}
bool EventSkimmer::checkMuonSelection(const edm::Event& iEvent) {
  edm::Handle<reco::MuonCollection> muons_;
  iEvent.getByLabel(_muonInputTag, muons_);
  bool decision = false;
  for ( reco::MuonCollection::const_iterator imuon = muons_->begin();
      imuon != muons_->end(); ++imuon ) {
    if ( imuon->pt() >= _muonPtCut)  {
      if (_verbosity > 0) std::cout << "Muon Selected with pt "<< imuon->pt() << " cutoff : " << _muonPtCut << std::endl;
      decision = true;
      break;
    }    
  }
  return decision;
}
bool EventSkimmer::checkElectronSelection(const edm::Event& iEvent) {
  edm::Handle<reco::ElectronCollection> electrons_;
  iEvent.getByLabel(_electronInputTag, electrons_);
  bool decision = false;
  for ( reco::ElectronCollection::const_iterator ielec = electrons_->begin();
      ielec != electrons_->end(); ++ielec ) {
    if ( ielec->pt() >= _electronPtCut)  {
      if (_verbosity > 0) std::cout << "Electron Selected with pt "<< ielec->pt() << " cutoff : " << _electronPtCut << std::endl;
      decision = true;
      break;
    }    
  }
  return decision;
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventSkimmer);
