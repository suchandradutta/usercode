#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TPRegexp.h"

#include "VHTauTau/TreeMaker/plugins/TriggerMuonBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

// Constructor
TriggerMuonBlock::TriggerMuonBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _hltInputTag(iConfig.getParameter<edm::InputTag>("hltInputTag")),
  _triggerEventTag(iConfig.getParameter<edm::InputTag>("triggerEventTag")),
  _muonInputTag(iConfig.getParameter<edm::InputTag>("muonSrcTag")),
  _muonMatchLabel(iConfig.getParameter<std::string>("muonMatchLabel")),
  _hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("hltPathsOfInterest")),
  _dMuonPath(iConfig.getParameter<std::string>("doubleMuonPath"))
{}
TriggerMuonBlock::~TriggerMuonBlock() {
}
void TriggerMuonBlock::beginJob() 
{
  std::string tree_name = "vhtree";
  TTree* tree = Utility::getTree(tree_name);
  cloneTriggerMuon = new TClonesArray("TriggerMuon");
  tree->Branch("TriggerMuon", &cloneTriggerMuon, 32000, 2);
  tree->Branch("nTriggerMuon", &fnTriggerMuon, "fnTriggerMuon/I");

  // Now book histograms
  edm::Service<TFileService> fileService;

  _histos1D["mass"] = fileService->make<TH1D>("mass",     "Mass_{Z} (GeV)",  90,   30., 120.);

  _histos1D["tagPt"]  = fileService->make<TH1D>("tagPt",  "p_{T} (GeV)",    100,    0., 100.);
  _histos1D["tagEta"] = fileService->make<TH1D>("tagEta", "#eta",  48,     -2.4,  2.4);
  _histos1D["tagPhi"] = fileService->make<TH1D>("tagPhi", "#phi",  60, -TMath::Pi(), TMath::Pi());

  _histos1D["tagTrigPt"]  = fileService->make<TH1D>("tagTrigPt",  "p_{T} (GeV)",    100,    0., 100.);
  _histos1D["tagTrigEta"] = fileService->make<TH1D>("tagTrigEta", "#eta",  48,     -2.4,  2.4);
  _histos1D["tagTrigPhi"] = fileService->make<TH1D>("tagTrigPhi", "#phi",  60, -TMath::Pi(), TMath::Pi());

  _histos1D["probePt"]  = fileService->make<TH1D>("probePt",  "p_{T} (GeV)",    100,    0., 100.);
  _histos1D["probeEta"] = fileService->make<TH1D>("probeEta", "#eta",  48,     -2.4,  2.4);
  _histos1D["probePhi"] = fileService->make<TH1D>("probePhi",   "#phi",  60, -TMath::Pi(), TMath::Pi());

  _histos1D["testPt"]   = fileService->make<TH1D>("testPt",   "p_{T} (GeV)",    100,    0., 100.);
  _histos1D["testEta"]  = fileService->make<TH1D>("testEta",  "#eta",  48,     -2.4,  2.4);
  _histos1D["testPhi"]  = fileService->make<TH1D>("testPhi",   "#phi",  60, -TMath::Pi(), TMath::Pi());

  // pt correlation plot
  _histos2D["ptTagCorr"] = fileService->make<TH2D>("ptTagCorr", "Object vs. candidate p_{T} (GeV)", 60, 0., 300., 60, 0., 300.);
  _histos2D["ptTagCorr"]->SetXTitle("candidate p_{T} (GeV)");
  _histos2D["ptTagCorr"]->SetYTitle("object p_{T} (GeV)");

  // eta correlation plot
  _histos2D["etaTagCorr"] = fileService->make<TH2D>("etaTagCorr", "Object vs. candidate #eta", 50, -2.5, 2.5, 50, -2.5, 2.5);
  _histos2D["etaTagCorr"]->SetXTitle("candidate #eta");
  _histos2D["etaTagCorr"]->SetYTitle("object #eta");

  // phi correlation plot
  _histos2D["phiTagCorr"] 
     = fileService->make<TH2D>("phiTagCorr", "Object vs. candidate #phi", 60, -TMath::Pi(), TMath::Pi(), 60, -TMath::Pi(), TMath::Pi());
  _histos2D["phiTagCorr"]->SetXTitle("candidate #phi");
  _histos2D["phiTagCorr"]->SetYTitle("object #phi");
}
void TriggerMuonBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed = true;
  if (hltConfig.init(iRun, iSetup, _hltInputTag.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("TriggerMuonBlock") << "HLT config with process name " 
                                     << _hltInputTag.process() << " successfully extracted";
  } 
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("TriggerMuonBlock") << "Error! HLT config extraction with process name " 
                                      << _hltInputTag.process() << " failed";
    // In this case, all access methods will return empty values!
  }
  const std::vector<std::string>& pathList = hltConfig.triggerNames();
  for (std::vector<std::string>::const_iterator it  = pathList.begin();
                                                it != pathList.end(); ++it) {
    std::string path = (*it);
    if (path.find(_dMuonPath) != std::string::npos) {
      _testPath = path;
      break;
    }
  }
}
void TriggerMuonBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneTriggerMuon->Clear();
  fnTriggerMuon = 0;  

  // trigger event
  edm::Handle<pat::TriggerEvent> triggerEvent;
  iEvent.getByLabel(_triggerEventTag, triggerEvent);

  // pat candidate collection
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(_muonInputTag, muons);

  // pat trigger helper to recieve for trigger
  // matching information
  const pat::helper::TriggerMatchHelper matchHelper;

  // loop over muon references for the tag muon
  bool matchDone = false;
  for ( size_t idxTag = 0; idxTag < muons->size(); ++idxTag) {
    if (matchDone) break;
    // single muon trigger match performed 
    // do we need to apply a loose set of offline cuts?
    const pat::TriggerObjectRef 
      trigRefTag(matchHelper.triggerMatchObject(muons, idxTag, _muonMatchLabel, iEvent, *triggerEvent));
    if (trigRefTag.isAvailable()) {
      double pt  = muons->at(idxTag).pt();
      double eta = muons->at(idxTag).eta();
      double phi = muons->at(idxTag).phi();
     
      _histos1D["tagPt"]->Fill(pt);
      _histos1D["tagEta"]->Fill(eta);
      _histos1D["tagPhi"]->Fill(phi);

      if (trigRefTag.isNonnull()) {             // check references (necessary!)
        double trigPt = trigRefTag->pt();
        double trigEta = trigRefTag->eta();
        double trigPhi = trigRefTag->phi();

        _histos1D["tagTrigPt"]->Fill(trigPt);
        _histos1D["tagTrigEta"]->Fill(trigEta);
        _histos1D["tagTrigPhi"]->Fill(trigPhi);

	_histos2D["ptTagCorr"]->Fill(pt, trigPt);
	_histos2D["etaTagCorr"]->Fill(eta, trigEta);
	_histos2D["phiTagCorr"]->Fill(phi, trigPhi);
      }

      // loop over muon references for the probe/test muon
      for (size_t idxProbe = 0; idxProbe < muons->size() && idxProbe != idxTag; ++idxProbe) {
        double mass = (muons->at(idxTag).p4() + muons->at(idxProbe).p4()).mass(); 
	_histos1D["mass"]->Fill(mass);
	if (fabs(mass-91.2) < 5) {
          double eta = muons->at(idxProbe).eta();
          double pt = muons->at(idxProbe).pt();
          double phi = muons->at(idxProbe).phi();
	  
	  _histos1D["probePt"]->Fill(pt);
	  _histos1D["probeEta"]->Fill(eta);
	  _histos1D["probePhi"]->Fill(phi);
	  
	  if (triggerEvent->path(_testPath)->wasRun() && triggerEvent->path(_testPath)->wasAccept()) {
	    _histos1D["testPt"]->Fill(pt);
	    _histos1D["testEta"]->Fill(eta);
	    _histos1D["testPhi"]->Fill(phi);
            matchDone = true;
          }   
	}
      }
    }
  }
  // get the trigger objects corresponding to the used matching (HLT muons) and
  // loop over selected trigger objects
  std::vector<pat::TriggerObjectRef> uniqueObjects;
  const pat::TriggerPathRefVector trigAllPaths(triggerEvent->pathRefs());
  for (pat::TriggerPathRefVector::const_iterator iPath  = trigAllPaths.begin();
       iPath != trigAllPaths.end(); iPath++) {
    std::string name = (**iPath).name();
    int nmatch = 0;
    for (std::vector<std::string>::const_iterator kt  = _hltPathsOfInterest.begin();
	 kt != _hltPathsOfInterest.end(); ++kt) {
      nmatch += TPRegexp(*kt).Match(name);
    }
    if (!nmatch) continue; 
    pat::TriggerObjectRefVector myObjects(triggerEvent->pathObjects( name, true));
    for (pat::TriggerObjectRefVector::const_iterator it  = myObjects.begin();
	 it != myObjects.end();	   ++it) {
      std::vector<pat::TriggerObjectRef>::iterator ifind = find(uniqueObjects.begin(), uniqueObjects.end(), (*it));
      if (ifind == uniqueObjects.end()) uniqueObjects.push_back((*it));
    }
  }
  if (_verbosity)  std::cout << " # of unique Trigger Objects " << uniqueObjects.size() << std::endl;

  for (std::vector<pat::TriggerObjectRef>::const_iterator it  = uniqueObjects.begin(); 
                                                   it != uniqueObjects.end(); 
                                                  ++it) 
  {
    if (fnTriggerMuon == kMaxTriggerMuon) {
      edm::LogInfo("TriggerMuonBlock") 
	<< "Too many Trigger Muons (HLT), fnTriggerMuon = " << fnTriggerMuon; 
      break;
    }
    triggerMuonB = new ((*cloneTriggerMuon)[fnTriggerMuon++]) TriggerMuon();
    triggerMuonB->eta = (**it).eta();
    triggerMuonB->phi = (**it).phi();
    triggerMuonB->pt  = (**it).pt();
    triggerMuonB->energy  = (**it).energy();
    
    pat::TriggerPathRefVector objectPaths(triggerEvent->objectPaths((*it), true));
    for (pat::TriggerPathRefVector::const_iterator iPath  = objectPaths.begin();
	                                             iPath != objectPaths.end();
                                                      ++iPath) 
    {
      std::string name = (**iPath).name();
      int nmatch = 0;
      for (std::vector<std::string>::const_iterator kt  = _hltPathsOfInterest.begin();
	   kt != _hltPathsOfInterest.end(); ++kt) {
	nmatch += TPRegexp(*kt).Match(name);
      }
      if (!nmatch) continue;
      triggerMuonB->pathList.push_back(name);
    }
    if (_verbosity) {
      std::cout << triggerMuonB->eta 
		<< ":" << triggerMuonB->phi 
		<< ":" << triggerMuonB->pt 
		<< ":" << triggerMuonB->energy
		<< std::endl;
      for (std::vector<std::string>::const_iterator jt  = triggerMuonB->pathList.begin();
	                                            jt != triggerMuonB->pathList.end(); ++jt)
      {
	std::cout << (*jt) << ",";
      }
      std::cout << std::endl;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerMuonBlock);
