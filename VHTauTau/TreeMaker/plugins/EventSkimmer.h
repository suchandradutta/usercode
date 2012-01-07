#ifndef __VHTauTau_TreeMaker_EventSkimmer_hh
#define __VHTauTau_TreeMaker_EventSkimmer_hh

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class EventSkimmer : public edm::EDFilter 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  bool filter(edm::Event& iEvent, edm::EventSetup const& iSetup);
  virtual void endJob();

public:
  explicit EventSkimmer(const edm::ParameterSet& iConfig);
  virtual ~EventSkimmer() {}
  bool checkMuonSelection(const edm::Event& iEvent);
  bool checkElectronSelection(const edm::Event& iEvent);

private:
  int _verbosity;
  int _totalEvents;
  int _selectedEvents;

  edm::InputTag _muonInputTag;
  edm::InputTag _electronInputTag;
  bool _doMuonSelection;
  bool _doElectronSelection;
  double _muonPtCut;
  double _electronPtCut;
  bool _createHisto;

};
#endif
