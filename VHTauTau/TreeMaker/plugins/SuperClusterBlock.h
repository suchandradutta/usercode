#ifndef __VHTauTau_TreeMaker_SuperClusterBlock_h
#define __VHTauTau_TreeMaker_SuperClusterBlock_h

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
class SuperCluster;

class SuperClusterBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit SuperClusterBlock(const edm::ParameterSet& iConfig);
  virtual ~SuperClusterBlock() {}

  enum {
    kMaxSuperCluster = 100
  };

private:
  TClonesArray* cloneSuperCluster; 
  int  fnSuperCluster;

  int _verbosity;

  edm::InputTag _ebInputTag;
  edm::InputTag _eeInputTag;
  edm::InputTag _ecalEBInputTag;
  edm::InputTag _ecalEEInputTag;
  edm::InputTag _trkInputTag;
  edm::InputTag _eleInputTag;

  vhtm::SuperCluster* scB;
};
#endif
