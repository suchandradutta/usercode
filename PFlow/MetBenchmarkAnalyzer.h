#ifndef __DQMOffline_PFTau_MetBenchmarkAnalyzer__
#define __DQMOffline_PFTau_MetBenchmarkAnalyzer__

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMOffline/PFTau/interface/MetBenchmark.h"


class MetBenchmarkAnalyzer: public edm::EDAnalyzer {
 public:
  
  MetBenchmarkAnalyzer(const edm::ParameterSet& parameterSet);
  
 private:
  void analyze(edm::Event const&, edm::EventSetup const&);
  void beginJob() ;
  void endJob();

  edm::InputTag matchLabel_;
  edm::InputTag inputLabel_;
  std::string benchmarkLabel_;
  
  MetBenchmark metBenchmark_;

};

#endif 
