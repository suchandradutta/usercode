#ifndef __VHTauTau_TreeMaker_MuonBlock_h
#define __VHTauTau_TreeMaker_MuonBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"

class TClonesArray;
class Muon;
class MuonMVAEstimator;

class MuonBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob(){}

public:
  explicit MuonBlock(const edm::ParameterSet& iConfig);
  virtual ~MuonBlock();
  void fillIsoDeposit(const pat::Muon& muo, vhtm::Muon* muonB);

  enum {
    kMaxMuon = 100
  };

private:
  TClonesArray* cloneMuon; 
  int  fnMuon;

  int _verbosity;

  edm::InputTag _muonInputTag;
  edm::InputTag _vtxInputTag;
  edm::InputTag _beamSpotInputTag;
  edm::InputTag _rhoInputTag;
  edm::InputTag _pfInputTag;
  bool _beamSpotCorr;
  std::string _muonID;

  vhtm::Muon* muonB;

  MuonMVAEstimator* fMuonIsoMVA;
  MuonMVAEstimator* fMuonIdMVA;
  MuonMVAEstimator* fMuonIsoRingsRadMVA;
  MuonEffectiveArea::MuonEffectiveAreaTarget target_;
};
#endif
