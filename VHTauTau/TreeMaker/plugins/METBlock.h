#ifndef __VHTauTau_TreeMaker_METBlock_hh
#define __VHTauTau_TreeMaker_METBlock_hh

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"

class TClonesArray;
class MET;

class METBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit METBlock(const edm::ParameterSet& iConfig);
  virtual ~METBlock() {}

  enum {
    kMaxMET = 100
  };

private:
  TClonesArray* cloneMET; 
  int  fnMET;

  int _verbosity;
  edm::InputTag _inputTag;

  vhtm::MET* metB;
};
#endif
