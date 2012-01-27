#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "TTree.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/JetBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

PFJetIDSelectionFunctor pfjetIDLoose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
PFJetIDSelectionFunctor pfjetIDTight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
pat::strbitset retpf = pfjetIDLoose.getBitTemplate();

// Constructor
JetBlock::JetBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("jetSrc")),
  _jecUncPath(iConfig.getParameter<std::string>("jecUncertainty")),
  _applyResJEC (iConfig.getParameter<bool>     ("applyResidualJEC")),
  _resJEC (iConfig.getParameter<std::string>   ("residualJEC"))
{}
void JetBlock::beginJob() 
{
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  cloneJet = new TClonesArray("vhtm::Jet");
  tree->Branch("Jet", &cloneJet, 32000, 2);
  tree->Branch("nJet", &fnJet, "fnJet/I");
}
void JetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneJet->Clear();
  fnJet = 0;

  JetCorrectionUncertainty *jecUnc = 0;
  JetCorrectorParameters *ResJetCorPar = 0;
  FactorizedJetCorrector *JEC = 0;
  bool applyResJECLocal = _applyResJEC;
  if (_applyResJEC) {
    try {
      edm::FileInPath fipUnc(_jecUncPath);;
      jecUnc = new JetCorrectionUncertainty(fipUnc.fullPath());

      edm::FileInPath fipRes(_resJEC);
      ResJetCorPar = new JetCorrectorParameters(fipRes.fullPath());
      std::vector<JetCorrectorParameters> vParam;
      vParam.push_back(*ResJetCorPar);
      JEC = new FactorizedJetCorrector(vParam);
    }
    catch (std::exception& ex) {
      edm::LogInfo("JetBlock") << "The following exception occurred:" << std::endl << ex.what();
      applyResJECLocal = false; 
    } 
  }
  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByLabel(_inputTag, jets);

  if (jets.isValid()) {
    edm::LogInfo("JetBlock") << "Total # PAT Jets: " << jets->size();
    for(std::vector<pat::Jet>::const_iterator it = jets->begin(); 
                                             it != jets->end(); ++it) {
      if (fnJet == kMaxJet) {
	edm::LogInfo("JetBlock") << "Too many PAT Jets, fnJet = " << fnJet; 
	break;
      }
      retpf.set(false);
      int passjetLoose = (pfjetIDLoose(*it, retpf)) ? 1 : 0;

      retpf.set(false);
      int passjetTight = (pfjetIDTight(*it, retpf)) ? 1 : 0;

      double corr = 1.;
      if (applyResJECLocal && iEvent.isRealData() ) {
        JEC->setJetEta(it->eta());
        JEC->setJetPt(it->pt()); // here you put the L2L3 Corrected jet pt
        corr = JEC->getCorrection();
      }

      if (jecUnc) {
        jecUnc->setJetEta(it->eta());
        jecUnc->setJetPt(it->pt()*corr); // the uncertainty is a function of the corrected pt
      }
      jetB = new ((*cloneJet)[fnJet++]) vhtm::Jet();

      // fill in all the vectors
      jetB->eta        = it->eta();
      jetB->phi        = it->phi();
      jetB->pt         = it->pt()*corr;
      jetB->pt_raw     = it->correctedJet("Uncorrected").pt();
      jetB->energy     = it->energy()*corr;
      jetB->energy_raw = it->correctedJet("Uncorrected").energy();
      jetB->jecUnc     = (jecUnc) ? jecUnc->getUncertainty(true) : -1;
      jetB->resJEC     = corr;
      jetB->partonFlavour               = it->partonFlavour();

      jetB->chargedEmEnergyFraction     = it->chargedEmEnergyFraction();
      jetB->chargedHadronEnergyFraction = it->chargedHadronEnergyFraction();
      jetB->chargedMuEnergyFraction     = it->chargedMuEnergyFraction();
      jetB->electronEnergyFraction      = it->electronEnergy()/it->energy();
      jetB->muonEnergyFraction          = it->muonEnergyFraction();
      jetB->neutralEmEnergyFraction     = it->neutralEmEnergyFraction();
      jetB->neutralHadronEnergyFraction = it->neutralHadronEnergyFraction();
      jetB->photonEnergyFraction        = it->photonEnergyFraction();

      jetB->chargedHadronMultiplicity   = it->chargedHadronMultiplicity();
      jetB->chargedMultiplicity         = it->chargedMultiplicity();
      jetB->electronMultiplicity        = it->electronMultiplicity();
      jetB->muonMultiplicity            = it->muonMultiplicity();
      jetB->neutralHadronMultiplicity   = it->neutralHadronMultiplicity();
      jetB->neutralMultiplicity         = it->neutralMultiplicity();
      jetB->photonMultiplicity          = it->photonMultiplicity();

      jetB->nConstituents               = it->numberOfDaughters();
      jetB->trackCountingHighEffBTag    = it->bDiscriminator("trackCountingHighEffBJetTags");
      jetB->trackCountingHighPurBTag    = it->bDiscriminator("trackCountingHighPurBJetTags");
      jetB->simpleSecondaryVertexHighEffBTag = it->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      jetB->simpleSecondaryVertexHighPurBTag = it->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      jetB->jetProbabilityBTag               = it->bDiscriminator("jetProbabilityBJetTags");
      jetB->jetBProbabilityBTag              = it->bDiscriminator("jetBProbabilityBJetTags");
      jetB->passLooseID = passjetLoose;
      jetB->passTightID = passjetTight;
      if (_verbosity > 0) 
      std::cout << "JetBlock: trackCountingHighEffBJetTag = " << jetB->trackCountingHighEffBTag
                << ", trackCountingHighPurBJetTag = " << jetB->trackCountingHighPurBTag
                << ", jetProbabilityBTag = " << jetB->jetProbabilityBTag  
                << ", jetBProbabilityBTag = " << jetB->jetBProbabilityBTag  
                << std::endl;
    }
  } 
  else {
    edm::LogError("JetBlock") << "Error! Can't get the product " << _inputTag;
  }
  if (jecUnc) delete jecUnc;
  if (ResJetCorPar) delete ResJetCorPar;
  if (JEC) delete JEC;
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetBlock);
