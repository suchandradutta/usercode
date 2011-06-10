#ifndef TauBlock_hh
#define TauBlock_hh

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <vector>

class TClonesArray;
class Tau;

class TauBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit TauBlock(const edm::ParameterSet& iConfig);
  virtual ~TauBlock() {}

  enum {
    kMaxTau = 100
  };

private:
  TClonesArray* cloneTau; 
  int  fnTau;
  int _verbosity;
  edm::InputTag _inputTag;

  Tau* tauB;
};
#endif
