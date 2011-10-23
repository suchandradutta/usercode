#ifndef __TreeMaker_TriggerMuonBlock_hh
#define __TreeMaker_TriggerMuonBlock_hh

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
class TriggerMuon;

class TriggerMuonBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit TriggerMuonBlock(const edm::ParameterSet& iConfig);
  virtual ~TriggerMuonBlock();

  enum {
    kMaxTriggerMuon = 100
  };

private:
  TClonesArray* cloneTriggerMuon; 
  int  fnTriggerMuon;

  int _verbosity;
  const edm::InputTag _hltInputTag;
  edm::InputTag _triggerEventTag;
  edm::InputTag _muonInputTag;
  std::string  _muonMatchLabel;
  const std::vector<std::string> _hltPathsOfInterest;
  const std::string _dMuonPath;

  TriggerMuon* triggerMuonB;

  std::string _testPath;

  // management of 1d histograms
  std::map< std::string, TH1D* > _histos1D;
  std::map< std::string, TH2D* > _histos2D;

  HLTConfigProvider hltConfig;
};
#endif
