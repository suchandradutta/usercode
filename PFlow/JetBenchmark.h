#ifndef DQMOffline_PFTau_JetBenchmark_h
#define DQMOffline_PFTau_JetBenchmark_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMOffline/PFTau/interface/Benchmark.h"
#include "DQMOffline/PFTau/interface/CandidateBenchmark.h"
#include "DQMOffline/PFTau/interface/MatchCandidateBenchmark.h"

#include "DataFormats/JetReco/interface/BasicJetCollection.h"

#include <vector>
class JetBenchmark : public Benchmark {

 public:

  JetBenchmark( float dRMax = 0.3,
		      bool matchCharge = true, 
		      Benchmark::Mode mode=Benchmark::DEFAULT) 
    : 
    Benchmark(mode), 
    candBench_(mode), matchCandBench_(mode), 
    dRMax_(dRMax), matchCharge_(matchCharge) {}
  
  virtual ~JetBenchmark();
  
  /// set the parameters
  void setParameters( const edm::ParameterSet& parameterSet);
  
  /// set directory (to use in ROOT)
  void setDirectory(TDirectory* dir);

  /// book histograms
  void setup();
  
  /// fill histograms with all particle
  template< class T, class C>
  void fill(const T& jetCollection,
	    const C& matchedJetCollection );


  void fillOne(const reco::Jet& jet,
	       const reco::Jet& matchedJet);

 protected:
  CandidateBenchmark      candBench_;
  MatchCandidateBenchmark matchCandBench_;

  TH2F*  delta_frac_VS_frac_muon_;
  TH2F*  delta_frac_VS_frac_photon_;
  TH2F*  delta_frac_VS_frac_electron_;
  TH2F*  delta_frac_VS_frac_charged_hadron_;
  TH2F*  delta_frac_VS_frac_neutral_hadron_;

  float dRMax_;
  bool  matchCharge_;
  bool  createPFractionHistos_;

};

#include "DQMOffline/PFTau/interface/Matchers.h"
template< class T, class C>
void JetBenchmark::fill(const T& jetCollection,
			const C& matchedJetCollection) {
  

  std::vector<int> matchIndices;
  PFB::match( jetCollection, matchedJetCollection, matchIndices, 
	      matchCharge_, dRMax_ );

  for (unsigned int i = 0; i < (jetCollection).size(); i++) {
    const reco::Jet& jet = jetCollection[i];

    if( !isInRange(jet.pt(), jet.eta(), jet.phi() ) ) continue;
    
    int iMatch = matchIndices[i];

    assert(iMatch< static_cast<int>(matchedJetCollection.size()));
 
    if( iMatch!=-1 ) {
      candBench_.fillOne(jet);
      matchCandBench_.fillOne(jet, matchedJetCollection[ iMatch ]);
      if (createPFractionHistos_) fillOne(jet, matchedJetCollection[ iMatch ]);
    }
  }
}
#endif 
