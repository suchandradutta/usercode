#include "DQMOffline/PFTau/plugins/JetBenchmarkAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//
// -- Constructor
//
JetBenchmarkAnalyzer::JetBenchmarkAnalyzer(const edm::ParameterSet& parameterSet)  
  
{
  inputLabel_          = parameterSet.getParameter<edm::InputTag>("InputCollection");
  matchLabel_          = parameterSet.getParameter<edm::InputTag>("MatchCollection");
  benchmarkLabel_      = parameterSet.getParameter<std::string>("BenchmarkLabel"); 

  jetBenchmark_.setParameters(parameterSet);  

}
//
// -- BeginJob
//
void JetBenchmarkAnalyzer::beginJob() {

  Benchmark::DQM_ = edm::Service<DQMStore>().operator->();
  // part of the following could be put in the base class
  std::string path = "ParticleFlow/" + benchmarkLabel_;
  Benchmark::DQM_->setCurrentFolder(path.c_str());
  std::cout<<"Histogram Folder path set to "<< path <<std::endl;
  jetBenchmark_.setup();  

}
//
// -- Analyze
//
void JetBenchmarkAnalyzer::analyze(edm::Event const& iEvent, 
				      edm::EventSetup const& iSetup) {
  edm::Handle< edm::View<reco::Jet> > jetCollection;
  iEvent.getByLabel(inputLabel_, jetCollection);   
  
  edm::Handle< edm::View<reco::Jet> > matchedJetCollection; 
  iEvent.getByLabel( matchLabel_, matchedJetCollection);

  if (jetCollection.isValid() && matchedJetCollection.isValid()) {
    jetBenchmark_.fill( *jetCollection, *matchedJetCollection );
  }
}

//
// -- EndJob
// 
void JetBenchmarkAnalyzer::endJob() {
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE (JetBenchmarkAnalyzer) ;
