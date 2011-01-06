#include "DQMOffline/PFTau/plugins/MetBenchmarkAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/METReco/interface/MET.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//
// -- Constructor
//
MetBenchmarkAnalyzer::MetBenchmarkAnalyzer(const edm::ParameterSet& parameterSet)  
  
{
  inputLabel_          = parameterSet.getParameter<edm::InputTag>("InputCollection");
  matchLabel_          = parameterSet.getParameter<edm::InputTag>("MatchCollection");
  benchmarkLabel_      = parameterSet.getParameter<std::string>("BenchmarkLabel"); 

  metBenchmark_.setParameters(parameterSet);  

}
//
// -- BeginJob
//
void MetBenchmarkAnalyzer::beginJob() {

  Benchmark::DQM_ = edm::Service<DQMStore>().operator->();
  // part of the following could be put in the base class
  std::string path = "ParticleFlow/" + benchmarkLabel_;
  Benchmark::DQM_->setCurrentFolder(path.c_str());
  std::cout<<"Histogram Folder path set to "<< path <<std::endl;
  metBenchmark_.setup();  

}
//
// -- Analyze
//
void MetBenchmarkAnalyzer::analyze(edm::Event const& iEvent, 
				      edm::EventSetup const& iSetup) {
  edm::Handle< edm::View<reco::MET> > metCollection;
  iEvent.getByLabel(inputLabel_, metCollection);   
  
  edm::Handle< edm::View<reco::MET> > matchedMetCollection; 
  iEvent.getByLabel( matchLabel_, matchedMetCollection);

  if (metCollection.isValid() && matchedMetCollection.isValid()) {
    metBenchmark_.fillOne( (*metCollection)[0], (*matchedMetCollection)[0] );
  }
}

//
// -- EndJob
// 
void MetBenchmarkAnalyzer::endJob() {
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (MetBenchmarkAnalyzer) ;
