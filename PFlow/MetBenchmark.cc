#include "DataFormats/METReco/interface/MET.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DQMOffline/PFTau/interface/Matchers.h"

#include "DQMOffline/PFTau/interface/MetBenchmark.h"

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

MetBenchmark::~MetBenchmark() {}


void MetBenchmark::setParameters( const edm::ParameterSet & parameterSet) {

  mode_           = (Benchmark::Mode) parameterSet.getParameter<int>( "mode" );
  variableEtBins_ = parameterSet.getParameter< std::vector<double> >( "VariableEtBins" );
  createMETSpecificHistos_ = parameterSet.getParameter<bool>( "CreateMETSpecificHistos" );
  setRange( parameterSet.getParameter<double>("ptMin"),
	    parameterSet.getParameter<double>("ptMax"),
	    parameterSet.getParameter<double>("etaMin"),
	    parameterSet.getParameter<double>("etaMax"),
	    parameterSet.getParameter<double>("phiMin"),
	    parameterSet.getParameter<double>("phiMax") );


  candBench_.setParameters(mode_);
  matchCandBench_.setParameters(mode_);
}

void MetBenchmark::setup() {
  candBench_.setup();
  matchCandBench_.setup();
  
  if (createMETSpecificHistos_) {
    float* etBins = new float[variableEtBins_.size()];
    for (size_t i = 0; i < variableEtBins_.size(); i++) {
      etBins[i] = variableEtBins_[i];
    }
    
    PhaseSpace dptPS      = PhaseSpace( 200, -500, 500);
    PhaseSpace setPS      = PhaseSpace( 50, 0.0, 3000);
    PhaseSpace dsetPS     = PhaseSpace( 50, -1000.0, 1000);
    PhaseSpace setOvsetPS = PhaseSpace( 100,0., 2.);
    
    delta_ex_VS_et_ = book2D("delta_ex_", 
			     "#DeltaME_{X}",variableEtBins_.size()-1, etBins,
			     dptPS.n, dptPS.m, dptPS.M );
    delta_ex_VS_et_ = book2D("delta_ey_", 
			     "#DeltaME_{Y}",variableEtBins_.size()-1, etBins,
			   dptPS.n, dptPS.m, dptPS.M );
    
    delta_set_VS_set_ = book2D("delta_set_VS_set_", 
			       ";SE_{T, true} (GeV);#DeltaSE_{T}",
			       setPS.n, setPS.m, setPS.M,
			       dsetPS.n, dsetPS.m, dsetPS.M );
    
    delta_set_Over_set_VS_set_ = book2D("delta_set_Over_set_VS_set_", 
					";SE_{T, true} (GeV);#DeltaSE_{T}/SE_{T}",
					setPS.n, setPS.m, setPS.M,
					setOvsetPS.n, setOvsetPS.m, setOvsetPS.M );
  }
}
void MetBenchmark::setDirectory(TDirectory* dir) {
  Benchmark::setDirectory(dir);

  candBench_.setDirectory(dir);
  matchCandBench_.setDirectory(dir);
}

void MetBenchmark::fillOne(const reco::MET& met,
			const reco::MET& matchedMet) {
  candBench_.fillOne(met);
  matchCandBench_.fillOne(met, matchedMet);
  /*  if (createMETSpecificHistos_) {
    if( !isInRange(met.pt(), met.eta(), met.phi() ) ) return;
    delta_ex_VS_et_->Fill(matchedMet.px(), met.px()-matchedMet.px());
    delta_ey_VS_et_->Fill(matchedMet.py(), met.py()-matchedMet.py());
    delta_set_VS_set_->Fill(matchedMet.sumEt(),met.sumEt()-matchedMet.sumEt());
    if ( matchedMet.sumEt()>0.001 )delta_set_Over_set_VS_set_->Fill(matchedMet.sumEt(),(met.sumEt()-matchedMet.sumEt())/matchedMet.sumEt());
    }*/
}
