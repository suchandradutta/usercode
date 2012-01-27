#ifndef __VHTauTau_TreeMaker_CommonVertexBlock_h
#define __VHTauTau_TreeMaker_CommonVertexBlock_h

#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"

class TransientTrackBuilder;
class TClonesArray;
class CommonVertex;

class CommonVertexBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit CommonVertexBlock(const edm::ParameterSet& iConfig);
  virtual ~CommonVertexBlock() {}

  void FitVertex(const TransientTrackBuilder* builder, vhtm::CommonVertex* v,
	   int nMuon, int nElectron, int nTau, std::string label);
private:
  TClonesArray* cloneCommonVertex; 

  int _verbosity;
  edm::InputTag _patMuonSrc;
  edm::InputTag _patElectronSrc;
  edm::InputTag _patTauSrc;

  double _minPtMuon;
  double _minPtElectron;
  double _minPtTau;

  std::vector<const pat::Muon*> _selectedMuons;
  std::vector<const pat::Electron*> _selectedElectrons;
  std::vector<const pat::Tau*> _selectedTaus;

  vhtm::CommonVertex* vertexB;
};
#endif