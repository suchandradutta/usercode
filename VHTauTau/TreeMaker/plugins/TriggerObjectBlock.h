#ifndef __TreeMaker_TriggerObjectBlock_hh
#define __TreeMaker_TriggerObjectBlock_hh

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
#include <map>

class TClonesArray;
class TriggerObject;

class TriggerObjectBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit TriggerObjectBlock(const edm::ParameterSet& iConfig);
  virtual ~TriggerObjectBlock();

  enum {
    kMaxTriggerObject = 100
  };

private:
  TClonesArray* cloneTriggerObject; 
  int  fnTriggerObject;

  int _verbosity;
  const edm::InputTag _hltInputTag;
  edm::InputTag _triggerEventTag;
  const std::string  _tagPathLabel;
  const std::string  _probePathLabel;
  const std::vector<std::string> _hltPathsOfInterest;

  TriggerObject* _triggerObject;

  HLTConfigProvider hltConfig;
};
#endif