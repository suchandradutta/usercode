#ifndef __TreeMaker_GenParticleBlock_hh
#define __TreeMaker_GenParticleBlock_hh

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <vector>

class TClonesArray;
class GenParticle;

class GenParticleBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit GenParticleBlock(const edm::ParameterSet& iConfig);
  virtual ~GenParticleBlock() {}

  enum {
    kMaxGenParticle = 1000
  };

private:
  TClonesArray* cloneGenParticle; 
  int  fnGenParticle;
  int _verbosity;
  edm::InputTag _inputTag;

  GenParticle* genParticleB;
};
#endif
