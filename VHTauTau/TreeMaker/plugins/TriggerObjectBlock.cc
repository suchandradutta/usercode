#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "VHTauTau/TreeMaker/plugins/TriggerObjectBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "TMath.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TPRegexp.h"

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
  TTree* tree = vhtm::Utility::getTree(tree_name);
  cloneTriggerObject = new TClonesArray("vhtm::TriggerObject");
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
  const pat::TriggerPathRefVector trigAllPaths(triggerEvent->pathRefs());
  pat::TriggerObjectRefVector myObjects(triggerEvent->objectRefs());
  int nObjects = 0; 
  for (pat::TriggerObjectRefVector::const_iterator it  = myObjects.begin();                                                 
       it != myObjects.end();    ++it) {                                 
    pat::TriggerPathRefVector myPaths = triggerEvent->objectPaths((*it));
    std::map <std::string, unsigned int> pathInfoMap;

    for (pat::TriggerPathRefVector::const_iterator ipath = myPaths.begin();
	 ipath != myPaths.end();    ++ipath) {
      std::string name = (**ipath).name();
      for (std::vector<std::string>::const_iterator kt  = _hltPathsOfInterest.begin();
	   kt != _hltPathsOfInterest.end(); ++kt) {
	std::string path_int = (*kt);
        if (name.find(path_int) == std::string::npos) continue;
	unsigned int val = 0;
	if (triggerEvent->path(name)->wasRun() && triggerEvent->path(name)->wasAccept()) val = 1;
	pathInfoMap.insert(std::pair<std::string, unsigned int> (name, val));
      }
    }
    if (pathInfoMap.size() > 0)  {
      if (fnTriggerObject == kMaxTriggerObject) {
	edm::LogInfo("TriggerObjectBlock") 
	  << "Too many Trigger Muons (HLT), fnTriggerObject = " << fnTriggerObject; 
	break;
      }
      nObjects++;
      _triggerObject = new ((*cloneTriggerObject)[fnTriggerObject++]) TriggerObject();
      _triggerObject->eta = (**it).eta();
      _triggerObject->phi = (**it).phi();
      _triggerObject->pt  = (**it).pt();
      _triggerObject->energy  = (**it).energy();
      for (std::map <std::string, unsigned int>::iterator imap = pathInfoMap.begin(); imap != pathInfoMap.end(); imap++) {
	_triggerObject->pathList.insert(std::pair<std::string, unsigned int> (imap->first, imap->second));
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
  if (_verbosity) std::cout << " # of Trigger Objects " << nObjects << std::endl;
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerObjectBlock);
