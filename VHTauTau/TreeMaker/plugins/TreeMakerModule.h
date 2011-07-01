#ifndef __TreeMaker_TreeMakerModule_hh
#define __TreeMaker_TreeMakerModule_hh

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <vector>

class TreeMakerModule : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob();

public:
  explicit TreeMakerModule(const edm::ParameterSet& iConfig);
  virtual ~TreeMakerModule() {}

private:
   int _verbosity;
  bool _createTree;
};
#endif
