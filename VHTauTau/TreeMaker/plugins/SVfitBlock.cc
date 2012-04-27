#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TVector3.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisFwd.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisBaseFwd.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitEventHypothesisByIntegration.h"
typedef std::vector<NSVfitEventHypothesisByIntegration> NSVfitEventHypothesisByIntegrationCollection;

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesisBase.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitSingleParticleHypothesis.h"

#include <Math/VectorUtil.h>
#include "VHTauTau/TreeMaker/plugins/SVfitBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

SVfitBlock::SVfitBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _diTauPairsSrc(iConfig.getParameter<edm::InputTag>("diTauPairs")),
  _genParticlesSrc(iConfig.getParameter<edm::InputTag>("genParticles"))
{
  //now do what ever initialization is needed
}

SVfitBlock::~SVfitBlock() {}
void SVfitBlock::beginJob()
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  cloneSVDiTau = new TClonesArray("vhtm::SVDiTau");
  tree->Branch("SVDiTau", &cloneSVDiTau, 32000, 2);
  tree->Branch("nSVDiTau", &fnSVDiTau,  "fnSVDiTau/I");
}
void SVfitBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<NSVfitEventHypothesisByIntegrationCollection> tauPairsInt;
  iEvent.getByLabel(_diTauPairsSrc, tauPairsInt);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(_genParticlesSrc, genParticles);

  const reco::GenParticle* muon = 0;
  const reco::GenParticle* tauPos = 0;
  const reco::GenParticle* tauNeg = 0;
  for (unsigned int j = 0; j < genParticles->size(); ++j) {
    const reco::GenParticle& gen = (*genParticles)[j];
    if (gen.status() != 3) continue;

    int pdgid = abs(gen.pdgId());
    double chg = gen.charge();
    if (pdgid == 13) { assert(!muon); muon = &gen; }
    if (pdgid == 15 && chg < 0) { assert(!tauNeg); tauNeg = &gen; }
    if (pdgid == 15 && chg > 0) { assert(!tauPos); tauPos = &gen; }
  }
  if (!muon || !tauPos || !tauNeg) return;

  for (unsigned i = 0 ; i < tauPairsInt->size() ; i++) {
    const NSVfitEventHypothesisBase& tauPair = (*tauPairsInt)[i];

    const NSVfitResonanceHypothesisBase* svFitResonanceHypothesis = tauPair.NSVfitEventHypothesisBase::resonance("Higgs");
    if (!svFitResonanceHypothesis) continue;

    const NSVfitResonanceHypothesisBase* svFitResonanceHypothesisW = tauPair.NSVfitEventHypothesisBase::resonance("W");
    if (!svFitResonanceHypothesisW) continue;

    const NSVfitSingleParticleHypothesisBase* svFitLeg1 = svFitResonanceHypothesis->daughter("leg1");
    if (!svFitLeg1) continue; 

    const NSVfitSingleParticleHypothesisBase* svFitLeg2 = svFitResonanceHypothesis->daughter("leg2");
    if (!svFitLeg2) continue; 

    const NSVfitSingleParticleHypothesisBase* svFitMuon = svFitResonanceHypothesisW->daughter("chargedLepton");
    if (!svFitMuon) continue;

    if (fnSVDiTau == kMaxSVDiTau) {
       edm::LogInfo("SVfitBlock") << "Too many SV DiTau, fnSVDiTau = " 
                                  << fnSVDiTau; 
       break;
    }
    svDiTauB = new ((*cloneSVDiTau)[fnSVDiTau++]) vhtm::SVDiTau();

    const NSVfitSingleParticleHypothesisBase* svFitLegNeg = (svFitLeg1->particle()->charge() < 0 ? svFitLeg1 : svFitLeg2);
    const NSVfitSingleParticleHypothesisBase* svFitLegPos = (svFitLeg1->particle()->charge() > 0 ? svFitLeg1 : svFitLeg2);
    svDiTauB->mass = svFitResonanceHypothesis->mass();
    svDiTauB->isValidSolution = svFitResonanceHypothesis->isValidSolution();

    if (ROOT::Math::VectorUtil::DeltaR(svFitLegNeg->particle()->p4(), tauNeg->p4()) < 0.1 &&
        ROOT::Math::VectorUtil::DeltaR(svFitLegPos->particle()->p4(), tauPos->p4()) < 0.1 &&
        ROOT::Math::VectorUtil::DeltaR(svFitMuon->particle()->p4(), muon->p4()) < 0.1)
    {
       svDiTauB->isGenMatched = true;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SVfitBlock);
