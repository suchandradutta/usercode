#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/TriggerBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"

static const unsigned int NmaxL1AlgoBit = 128;
static const unsigned int NmaxL1TechBit = 64;

// Constructor
TriggerBlock::TriggerBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _l1InputTag(iConfig.getParameter<edm::InputTag>("L1InputTag")),
  _hltInputTag(iConfig.getParameter<edm::InputTag>("HLTInputTag")),
  _hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("HLTPathsOfInterest"))
{}
void TriggerBlock::beginJob() 
{
  std::string tree_name = "vhtree";
  TTree* tree = Utility::getTree(tree_name);
  cloneTrigger = new TClonesArray("Trigger");
  tree->Branch("Trigger", &cloneTrigger, 32000, 2);
}
void TriggerBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed = true;
  if (hltConfig.init(iRun, iSetup, _hltInputTag.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("TriggerBlock") << "HLT config with process name " 
                                 << _hltInputTag.process() << " successfully extracted";
  } 
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("TriggerBlock") << "Error! HLT config extraction with process name " 
                                  << _hltInputTag.process() << " failed";
    // In this case, all access methods will return empty values!
  }
}
void TriggerBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneTrigger->Clear();

  // Create Trigger Object
  triggerB = new ( (*cloneTrigger)[0] ) Trigger();

  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByLabel(_l1InputTag, l1GtReadoutRecord);

  if (l1GtReadoutRecord.isValid()) {
    edm::LogInfo("TriggerBlock") << "Successfully obtained " << _l1InputTag;
    for (unsigned int i = 0; i < NmaxL1AlgoBit; ++i) {
      triggerB->l1physbits.push_back(l1GtReadoutRecord->decisionWord()[i] ? 1 : 0);
    }
    for (unsigned int i = 0; i < NmaxL1TechBit; ++i) {
      triggerB->l1techbits.push_back( l1GtReadoutRecord->technicalTriggerWord()[i] ? 1 : 0 );
    }
  } 
  else {
    edm::LogError("TriggerBlock") << "Error! Can't get the product " << _l1InputTag;
  }
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(_hltInputTag, triggerResults);
  if (triggerResults.isValid()) {
    edm::LogInfo("TriggerBlock") << "Successfully obtained " << _hltInputTag;
    for (unsigned int i = 0; i < triggerResults->size(); i++) {
      triggerB->hltbits.push_back(triggerResults->at(i).accept() ? 1 : 0);
    }
    for (std::vector<std::string>::const_iterator it = _hltPathsOfInterest.begin();
                                                 it != _hltPathsOfInterest.end(); ++it) {
      int fired = 0;
      unsigned int index = hltConfig.triggerIndex(*it);
      if (index < triggerResults->size()) {
        if (triggerResults->accept(index)) fired = 1;
      } 
      else {
	edm::LogInfo("TriggerBlock") << "Requested HLT path \"" << (*it) << "\" does not exist";
      }
      triggerB->hltresults.push_back(fired) ;

      int prescale = -1;
      if (hltConfig.prescaleSet(iEvent, iSetup) < 0) {
	edm::LogError("TriggerBlock") << "Error! The prescale set index number could not be obtained";
      } 
      else {
        prescale = hltConfig.prescaleValue(iEvent, iSetup, *it);
      }
      triggerB->hltprescales.push_back(prescale);
    }
  } 
  else {
    edm::LogError("TriggerBlock") << "Error! Can't get the product " << _hltInputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerBlock);
