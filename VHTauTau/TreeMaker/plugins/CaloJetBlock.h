#ifndef CaloJetBlock_hh
#define CaloJetBlock_hh

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include <string>
#include <vector>

class TClonesArray;
class CaloJet;

class CaloJetBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit CaloJetBlock(const edm::ParameterSet& iConfig);
  virtual ~CaloJetBlock() {}

  enum {
    kMaxCaloJet = 100
  };
private:
  TClonesArray* cloneCaloJet; 
  int fnCaloJet;

  int _verbosity;
  edm::InputTag _inputTag;
  double _electronPt;
  double _electronIso;
  double _muonPt;
  double _muonIso;
  std::string _jecUncPath;
  bool _applyResJEC;
  std::string _resJEC;


  CaloJet* caloJetB;
};
#endif
