#ifndef __DQMOffline_PFTau_JetBenchmarkAnalyzer__
#define __DQMOffline_PFTau_JetBenchmarkAnalyzer__

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMOffline/PFTau/interface/JetBenchmark.h"


class JetBenchmarkAnalyzer: public edm::EDAnalyzer {
 public:
  
  JetBenchmarkAnalyzer(const edm::ParameterSet& parameterSet);
  
 private:
  void analyze(edm::Event const&, edm::EventSetup const&);
  void beginJob() ;
  void endJob();

  edm::InputTag matchLabel_;
  edm::InputTag inputLabel_;
  std::string benchmarkLabel_;
  
  JetBenchmark jetBenchmark_;

};

#endif 
