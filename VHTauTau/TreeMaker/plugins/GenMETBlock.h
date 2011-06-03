#ifndef GenMETBlock_hh
#define GenMETBlock_hh

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <vector>

class GenMET;
class TClonesArray;

class GenMETBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit GenMETBlock(const edm::ParameterSet& iConfig);
  virtual ~GenMETBlock() {}

  enum {
    kMaxGenMET = 100
  };

private:
  TClonesArray* cloneGenMET; 
  int  fnGenMET;
  int _verbosity;
  edm::InputTag _inputTag;

  GenMET* genMetB;
};
#endif
