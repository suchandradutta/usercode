#ifndef DQMOffline_PFTau_JetBenchmark_h
#define DQMOffline_PFTau_JetBenchmark_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMOffline/PFTau/interface/Benchmark.h"
#include "DQMOffline/PFTau/interface/CandidateBenchmark.h"
#include "DQMOffline/PFTau/interface/MatchCandidateBenchmark.h"
#include "DataFormats/METReco/interface/METCollection.h"

#include <vector>
class MetBenchmark : public Benchmark {

 public:

  MetBenchmark( Benchmark::Mode mode=Benchmark::DEFAULT) 
    : 
    Benchmark(mode), 
    candBench_(mode), matchCandBench_(mode) {} 
  
  virtual ~MetBenchmark();
  
  /// set the parameters
  void setParameters( const edm::ParameterSet& parameterSet);
  
  /// set directory (to use in ROOT)
  void setDirectory(TDirectory* dir);

  /// book histograms
  void setup();
  

  void fillOne(const reco::MET& met,
	       const reco::MET& matchedMet);

 protected:

  TH2F*   delta_ex_VS_et_;
  TH2F*   delta_ey_VS_et_;
  TH2F*   delta_set_VS_set_;
  TH2F*   delta_set_Over_set_VS_set_;

  std::vector<double> variableEtBins_;
 
  CandidateBenchmark      candBench_;
  MatchCandidateBenchmark matchCandBench_;
  bool  createMETSpecificHistos_;

};
#endif 
