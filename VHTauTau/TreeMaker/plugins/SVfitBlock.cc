#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TMath.h"

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
  // Reset the TClonesArray and the nObj variables
  cloneSVDiTau->Clear();
  fnSVDiTau = 0;

  edm::Handle<NSVfitEventHypothesisByIntegrationCollection> tauPairsInt;
  iEvent.getByLabel(_diTauPairsSrc, tauPairsInt);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(_genParticlesSrc, genParticles);

  if (_verbosity) std::cout << "nTauPair = " << tauPairsInt->size() << std::endl;

  for (NSVfitEventHypothesisByIntegrationCollection::const_iterator it  = tauPairsInt->begin(); 
                                                                    it != tauPairsInt->end(); ++it) {

    const NSVfitEventHypothesisBase& tauPair = (*it);

    const NSVfitResonanceHypothesisBase* svFitResonanceHypothesis 
      = tauPair.NSVfitEventHypothesisBase::resonance("Higgs");
    if (!svFitResonanceHypothesis) continue;

    const NSVfitResonanceHypothesisBase* svFitResonanceHypothesisW 
      = tauPair.NSVfitEventHypothesisBase::resonance("W");
    if (!svFitResonanceHypothesisW) continue;

    const NSVfitSingleParticleHypothesisBase* svFitLeg1 
      = svFitResonanceHypothesis->daughter("leg1");
    if (!svFitLeg1) continue; 

    const NSVfitSingleParticleHypothesisBase* svFitLeg2 
      = svFitResonanceHypothesis->daughter("leg2");
    if (!svFitLeg2) continue; 

    const NSVfitSingleParticleHypothesisBase* svFitMuon 
      = svFitResonanceHypothesisW->daughter("chargedLepton");
    if (!svFitMuon) continue;

    if (fnSVDiTau == kMaxSVDiTau) {
      edm::LogInfo("SVfitBlock") << "Too many SV DiTau, fnSVDiTau = " << fnSVDiTau; 
      break;
    }
    svDiTauB = new ((*cloneSVDiTau)[fnSVDiTau++]) vhtm::SVDiTau();

    const NSVfitSingleParticleHypothesisBase* svFitLegNeg 
      = (svFitLeg1->particle()->charge() < 0) ? svFitLeg1 : svFitLeg2;
    const NSVfitSingleParticleHypothesisBase* svFitLegPos 
      = (svFitLeg1->particle()->charge() > 0) ? svFitLeg1 : svFitLeg2;

    svDiTauB->LegNegPx = svFitLegNeg->particle()->px();
    svDiTauB->LegNegPy = svFitLegNeg->particle()->py();
    svDiTauB->LegNegPz = svFitLegNeg->particle()->pz();
    svDiTauB->LegNegE  = svFitLegNeg->particle()->energy();
    svDiTauB->LegPosPx = svFitLegPos->particle()->px();
    svDiTauB->LegPosPy = svFitLegPos->particle()->py();
    svDiTauB->LegPosPz = svFitLegPos->particle()->pz();
    svDiTauB->LegPosE  = svFitLegPos->particle()->energy();
    svDiTauB->MuonPx   = svFitMuon->particle()->px();
    svDiTauB->MuonPy   = svFitMuon->particle()->py();
    svDiTauB->MuonPz   = svFitMuon->particle()->pz();
    svDiTauB->MuonE    = svFitMuon->particle()->energy();
    svDiTauB->mass     = svFitResonanceHypothesis->mass();
    svDiTauB->isValidSolution = svFitResonanceHypothesis->isValidSolution();
    svDiTauB->sigmaUp = svFitResonanceHypothesis->massErrUp();
    svDiTauB->sigmaDn = svFitResonanceHypothesis->massErrDown();

    if (_verbosity)
    std::cout << "SVfitBlock mass:" << std::setw(9) << svDiTauB->mass 
              << " isValidSolution: " << svDiTauB->isValidSolution 
              << std::endl;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SVfitBlock);
