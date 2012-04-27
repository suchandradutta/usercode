
#ifndef __VHTauTau_TreeMaker_SVfitBlock_h
#define __VHTauTau_TreeMaker_SVfitBlock_h

#include <string>
#include <vector>

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"

class TClonesArray;
class SVDiTau;

class SVfitBlock : public edm::EDAnalyzer {
public:
  explicit SVfitBlock(const edm::ParameterSet&);
  ~SVfitBlock();

  enum {
    kMaxSVDiTau = 10
  };

private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() {}

  TClonesArray* cloneSVDiTau; 
  int fnSVDiTau;

  int _verbosity;
  edm::InputTag _diTauPairsSrc;
  edm::InputTag _genParticlesSrc;

  vhtm::SVDiTau* svDiTauB;
};
#endif
