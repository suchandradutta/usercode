#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TMath.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TPRegexp.h"

#include "VHTauTau/TreeMaker/plugins/TriggerObjectBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

// Constructor
TriggerObjectBlock::TriggerObjectBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _hltInputTag(iConfig.getParameter<edm::InputTag>("hltInputTag")),
  _triggerEventTag(iConfig.getParameter<edm::InputTag>("triggerEventTag")),
  _hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("hltPathsOfInterest"))
{}
TriggerObjectBlock::~TriggerObjectBlock() {
}
void TriggerObjectBlock::beginJob() 
{
  std::string tree_name = "vhtree";
  TTree* tree = Utility::getTree(tree_name);
  cloneTriggerObject = new TClonesArray("TriggerObject");
  tree->Branch("TriggerObject", &cloneTriggerObject, 32000, 2);
  tree->Branch("nTriggerObject", &fnTriggerObject, "fnTriggerObject/I");

  // Now book histograms
  edm::Service<TFileService> fileService;

}
void TriggerObjectBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed = true;
  if (hltConfig.init(iRun, iSetup, _hltInputTag.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("TriggerObjectBlock") << "HLT config with process name " 
                                     << _hltInputTag.process() << " successfully extracted";
  } 
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("TriggerObjectBlock") << "Error! HLT config extraction with process name " 
                                      << _hltInputTag.process() << " failed";
    // In this case, all access methods will return empty values!
  }
}
void TriggerObjectBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneTriggerObject->Clear();
  fnTriggerObject = 0;  

  // trigger event
  edm::Handle<pat::TriggerEvent> triggerEvent;
  iEvent.getByLabel(_triggerEventTag, triggerEvent);

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
    if (fnTriggerObject == kMaxTriggerObject) {
      edm::LogInfo("TriggerObjectBlock") 
	<< "Too many Trigger Muons (HLT), fnTriggerObject = " << fnTriggerObject; 
      break;
    }
    _triggerObject = new ((*cloneTriggerObject)[fnTriggerObject++]) TriggerObject();
    _triggerObject->eta = (**it).eta();
    _triggerObject->phi = (**it).phi();
    _triggerObject->pt  = (**it).pt();
    _triggerObject->energy  = (**it).energy();
    
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
      unsigned int val = 0;
      if (triggerEvent->path(name)->wasRun() && triggerEvent->path(name)->wasAccept()) val = 1;
      _triggerObject->pathList.insert(std::pair<std::string, unsigned int> (name, val));
    }
    if (_verbosity) {
      std::cout << _triggerObject->eta 
		<< ":" << _triggerObject->phi 
		<< ":" << _triggerObject->pt 
		<< ":" << _triggerObject->energy
		<< std::endl;
      for (std::map<std::string, unsigned int>::const_iterator jt  = _triggerObject->pathList.begin();
	                                            jt != _triggerObject->pathList.end(); ++jt)
      {
	std::cout << jt->first << " flag " << jt->second << std::endl;
      }
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerObjectBlock);
