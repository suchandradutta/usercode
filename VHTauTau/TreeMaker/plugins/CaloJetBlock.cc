#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "VHTauTau/TreeMaker/interface/Utility.h"
#include "VHTauTau/TreeMaker/plugins/CaloJetBlock.h"

JetIDSelectionFunctor jetIDLoose(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE);
JetIDSelectionFunctor jetIDTight(JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::TIGHT);

pat::strbitset ret = jetIDLoose.getBitTemplate();

// Constructor
CaloJetBlock::CaloJetBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("caloJetSrc")),
  _electronPt (iConfig.getParameter<double>    ("electronPt")),
  _electronIso (iConfig.getParameter<double>   ("electronIso")),
  _muonPt (iConfig.getParameter<double>        ("muonPt")),
  _muonIso (iConfig.getParameter<double>       ("muonIso")),
  _jecUncPath(iConfig.getParameter<std::string>("jecUncertainty")),
  _applyResJEC (iConfig.getParameter<bool>     ("applyResidualJEC")),
  _resJEC (iConfig.getParameter<std::string>   ("residualJEC"))
{}
void CaloJetBlock::beginJob() 
{
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  cloneCaloJet = new TClonesArray("vhtm::CaloJet");
  tree->Branch("CaloJet", &cloneCaloJet, 32000, 2);
  tree->Branch("nCaloJet", &fnCaloJet, "fnCaloJet/I");
}
void CaloJetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneCaloJet->Clear();
  fnCaloJet = 0;

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
      edm::LogInfo("CaloJetBlock") << "The following exception occurred:" 
                                   << std::endl 
                                   << ex.what();
      applyResJECLocal = false; 
    } 
  }

  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByLabel(_inputTag, jets);

  if (jets.isValid()) {
    edm::LogInfo("CaloJetBlock") << "Total # CaloJets: " << jets->size();
    for (std::vector<pat::Jet>::const_iterator it = jets->begin(); 
                                              it != jets->end(); ++it) {
      if (fnCaloJet == kMaxCaloJet) {
	edm::LogInfo("CaloJetBlock") << "Too many PAT Jets, fnCaloJet = " << fnCaloJet; 
	break;
      }

      ret.set(false);
      int passjetLoose = (jetIDLoose(*it, ret)) ? 1 : 0;

      ret.set(false);
      int passjetTight = (jetIDTight(*it, ret)) ? 1 : 0;

      int ovrlps = 0;
      /* overlaps with good electrons (with different electron IDs) and muons are handled bitwise
         bit 0: eidRobustLoose
         bit 1: eidRobustTight
         bit 2: eidLoose
         bit 3: eidTight
         bit 4: eidRobustHighEnergy
         bit 5: HEEPId
         bit 6: GlobalMuonPromptTight
      */
      const reco::CandidatePtrVector & electrons = it->overlaps("electrons");
      for (size_t i = 0; i < electrons.size(); ++i) {
        // try to cast into pat::Electron
        const pat::Electron *electron = dynamic_cast<const pat::Electron *>(&*electrons[i]);
        if (electron) {
          double tiso = electron->trackIso() + electron->ecalIso() + electron->hcalIso();
          double pt   = electron->pt();
          double ratio = (pt > 0) ? tiso/pt : 0;
	  if (electron->electronID("simpleEleId70cIso")>0. 
                && ratio < _electronIso
	        && pt >_electronPt) ovrlps |= 1<<0;
          if (electron->electronID("simpleEleId80cIso")>0. 
                && ratio < _electronIso
	        && pt >_electronPt) ovrlps |= 1<<1;
          if (electron->electronID("simpleEleId85cIso")>0.       
                && ratio < _electronIso
	        && pt >_electronPt) ovrlps |= 1<<2;
          if (electron->electronID("simpleEleId90cIso")>0.       
                && ratio <  _electronIso
	        && pt >_electronPt) ovrlps |= 1<<3;
          if (electron->electronID("simpleEleId95cIso")>0. 
                && ratio < _electronIso
	        && pt >_electronPt) ovrlps |= 1<<4;
          if (electron->userInt("HEEPId")==0
	        && pt > _electronPt) ovrlps |= 1<<5;
        }
      }
      const reco::CandidatePtrVector & muons = it->overlaps("muons");
      for (size_t i = 0; i < muons.size(); ++i) {
        // try to cast into pat::Muon
        const pat::Muon *muon = dynamic_cast<const pat::Muon *>(&*muons[i]);
        if (muon) {
          double tiso = muon->trackIso() + muon->ecalIso() + muon->hcalIso();
          double pt   = muon->pt();
          double ratio = (pt > 0) ? tiso/pt : 0.0;
	  if (muon->muonID("GlobalMuonPromptTight")
	      && ratio < _muonIso
	      && pt > _muonPt ) ovrlps |= 1<<6;
        }
      }
      double corr = 1.;
      if (applyResJECLocal && iEvent.isRealData()) {
        JEC->setJetEta(it->eta());
        JEC->setJetPt(it->pt()); // here you put the L2L3 Corrected jet pt
        corr = JEC->getCorrection();
      }

      if (jecUnc) {
        jecUnc->setJetEta(it->eta());
        jecUnc->setJetPt(it->pt()*corr); // the uncertainty is a function of the corrected pt
      }

      caloJetB = new ((*cloneCaloJet)[fnCaloJet++]) vhtm::CaloJet();

      // fill in all the vectors
      caloJetB->eta        = it->eta();
      caloJetB->phi        = it->phi();
      caloJetB->pt         = it->pt() * corr;
      caloJetB->pt_raw     = it->correctedJet("Uncorrected").pt();
      caloJetB->energy     = it->energy() * corr;
      caloJetB->energy_raw = it->correctedJet("Uncorrected").energy();
      caloJetB->jecUnc     = (jecUnc) ? jecUnc->getUncertainty(true) : -1;
      caloJetB->resJEC     = corr;
      caloJetB->overlaps   = ovrlps;
      caloJetB->partonFlavour = it->partonFlavour();
      caloJetB->emf           = it->emEnergyFraction();
      caloJetB->resEmf        = it->jetID().restrictedEMF;
      caloJetB->hadf          = it->energyFractionHadronic();
      caloJetB->n90Hits       = it->jetID().n90Hits;
      caloJetB->fHPD          = it->jetID().fHPD;
      caloJetB->fRBX          = it->jetID().fRBX;
      caloJetB->sigmaEta      = std::sqrt(it->etaetaMoment());
      caloJetB->sigmaPhi      = std::sqrt(it->phiphiMoment());
      caloJetB->trackCountingHighEffBTag         = it->bDiscriminator("trackCountingHighEffBJetTags");
      caloJetB->trackCountingHighPurBTag         = it->bDiscriminator("trackCountingHighPurBJetTags");
      caloJetB->simpleSecondaryVertexHighEffBTag = it->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
      caloJetB->simpleSecondaryVertexHighPurBTag = it->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
      caloJetB->jetProbabilityBTag               = it->bDiscriminator("jetProbabilityBJetTags");
      caloJetB->jetBProbabilityBTag              = it->bDiscriminator("jetBProbabilityBJetTags");
      caloJetB->combinedSecondaryVertexBTag      = it->bDiscriminator("combinedSecondaryVertexBJetTags");
      caloJetB->combinedSecondaryVertexMVABTag   = it->bDiscriminator("combinedSecondaryVertexMVABJetTags");
      caloJetB->combinedInclusiveSecondaryVertexBTag = it->bDiscriminator("combinedInclusiveSecondaryVertexBJetTags");
      caloJetB->combinedMVABTag                  = it->bDiscriminator("combinedMVABJetTags");
      caloJetB->passLooseID = passjetLoose;
      caloJetB->passTightID = passjetTight;
    }
  } 
  else {
    edm::LogError("CaloJetBlock") << "Error >> Failed to get pat::Jet collection for label = " 
                                  << _inputTag;
  }
  if (jecUnc) delete jecUnc;
  if (ResJetCorPar) delete ResJetCorPar;
  if (JEC) delete JEC;
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CaloJetBlock);
