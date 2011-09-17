#ifndef __TreeMaker_EventBlock_hh
#define __TreeMaker_EventBlock_hh

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <string>
#include <vector>

class Event;

class EventBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
  virtual void endJob(){}

public:
  explicit EventBlock(const edm::ParameterSet& iConfig);
  virtual ~EventBlock();

private:
  TClonesArray* cloneEvent; 
  int _verbosity;

  const edm::InputTag   _l1InputTag;
  const edm::InputTag   _vtxInputTag;
  const edm::InputTag   _trkInputTag;
  const edm::InputTag   _hcalNoiseInputTag;
  const unsigned int    _vtxMinNDOF;
  const double          _vtxMaxAbsZ, _vtxMaxd0;
  const unsigned int    _numTracks;
  const double          _hpTrackThreshold;

  std::vector<int> *_nPU;
  std::vector<int> *_bunchCrossing;

  Event* eventB;
};
#endif
