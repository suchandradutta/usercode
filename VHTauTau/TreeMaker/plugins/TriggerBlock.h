#ifndef TriggerBlock_hh
#define TriggerBlock_hh

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include <string>
#include <vector>

class TTree;
class TClonesArray;
class Trigger;

class TriggerBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit TriggerBlock(const edm::ParameterSet& iConfig);
  virtual ~TriggerBlock() {}

private:
  TClonesArray* cloneTrigger; 

  TTree* _tree;
  int _verbosity;
  edm::InputTag _inputTag;

  const edm::InputTag _l1InputTag;
  const edm::InputTag _hltInputTag;
  const std::vector<std::string> _hltPathsOfInterest;
  HLTConfigProvider hltConfig;

  Trigger* triggerB;
};
#endif
