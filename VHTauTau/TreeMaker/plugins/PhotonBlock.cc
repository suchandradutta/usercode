#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "TClonesArray.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

#include "VHTauTau/TreeMaker/plugins/PhotonBlock.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

// Constructor
PhotonBlock::PhotonBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _treeName(iConfig.getParameter<std::string>("treeName")),
  _photonInputTag(iConfig.getParameter<edm::InputTag>("photonSrc"))
{}
void PhotonBlock::beginJob() 
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree(_treeName);
  clonePhoton = new TClonesArray("vhtm::Photon");
  tree->Branch("Photon", &clonePhoton, 32000, 2);
  tree->Branch("nPhoton", &fnPhoton,  "fnPhoton/I");
}
void PhotonBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  clonePhoton->Clear();
  fnPhoton = 0;

  edm::Handle<std::vector<pat::Photon> > photons;
  iEvent.getByLabel(_photonInputTag, photons);

  if (photons.isValid()) {
    edm::LogInfo("PhotonBlock") << "Total # PAT Photons: " << photons->size();
    for (std::vector<pat::Photon>::const_iterator it  = photons->begin(); 
                                                  it != photons->end(); 
                                                ++it) {
      if (fnPhoton == kMaxPhoton) {
	edm::LogInfo("PhotonBlock") << "Too many PAT Photon, fnPhoton = " 
                                    << fnPhoton; 
	break;
      }
      photonB = new ((*clonePhoton)[fnPhoton++]) vhtm::Photon();
      photonB->et     = it->et();
      photonB->eta    = it->eta();
      photonB->phi    = it->phi();
      photonB->energy = it->energy();
      photonB->theta  = it->theta();
      photonB->vx     = it->vx();
      photonB->vy     = it->vy();
      photonB->vz     = it->vz();

      const reco::SuperClusterRef sCluster = it->superCluster(); 
      photonB->scEnergy    = sCluster->energy();
      photonB->scEta       = sCluster->eta();
      photonB->scPhi       = sCluster->phi();
      photonB->scSize      = sCluster->clustersSize();
      photonB->scEtaWidth  = sCluster->etaWidth();
      photonB->scPhiWidth  = sCluster->etaWidth();
      photonB->scEt        = sCluster->energy()/cosh(sCluster->eta());
      photonB->scRawEnergy = sCluster->rawEnergy();
      photonB->scx         = sCluster->x();
      photonB->scx         = sCluster->y();
      photonB->scx         = sCluster->z();

      photonB->isoEcalRecHit03    = it->ecalRecHitSumEtConeDR03();
      photonB->isoHcalRecHit03    = it->hcalTowerSumEtConeDR03();
      photonB->isoSolidTrkCone03  = it->trkSumPtSolidConeDR03();
      photonB->isoHollowTrkCone03 = it->trkSumPtHollowConeDR03();
      photonB->nTrkSolidCone03    = it->nTrkSolidConeDR03();
      photonB->nTrkHollowCone03   = it->nTrkHollowConeDR03();

      photonB->isoEcalRecHit04    = it->ecalRecHitSumEtConeDR04();
      photonB->isoHcalRecHit04    = it->hcalTowerSumEtConeDR04();
      photonB->isoSolidTrkCone04  = it->trkSumPtSolidConeDR04();
      photonB->isoHollowTrkCone04 = it->trkSumPtHollowConeDR04();
      photonB->nTrkSolidCone04    = it->nTrkSolidConeDR04();
      photonB->nTrkHollowCone04   = it->nTrkHollowConeDR04();

      photonB->hasPixelSeed       = it->hasPixelSeed(); 
      photonB->ecalIso            = it->ecalIso();
      photonB->hcalIso            = it->hcalIso();
      photonB->trackIso           = it->trackIso();

      photonB->chargedHadIso      = it->chargedHadronIso();
      photonB->neutralHadIso      = it->neutralHadronIso();
      photonB->photonIso          = it->photonIso();


      int fidFlag = 0;
      if (it->isEB())        fidFlag |= (1 << 0);
      if (it->isEE())        fidFlag |= (1 << 1);
      if (it->isEBEtaGap())  fidFlag |= (1 << 2);
      if (it->isEBPhiGap())  fidFlag |= (1 << 3);
      if (it->isEERingGap()) fidFlag |= (1 << 4);
      if (it->isEEDeeGap())  fidFlag |= (1 << 5);
      if (it->isEBEEGap())   fidFlag |= (1 << 6);
      photonB->fidFlag = fidFlag;

      photonB->isEB               = it->isEB() ? true : false;
      photonB->isEE               = it->isEE() ? true : false;
      photonB->isEBGap            = it->isEBGap() ? true : false ;
      photonB->isEEGap            = it->isEEGap() ? true : false;
      photonB->isEBEEGap          = it->isEBEEGap() ? true : false;

      photonB->r9                 = it->r9();
      photonB->hoe                = it->hadronicOverEm();
      photonB->sigmaEtaEta        = it->sigmaEtaEta();
      photonB->sigmaIEtaIEta      = it->sigmaIetaIeta();
      photonB->e1x5               = it->e1x5();
      photonB->e2x5               = it->e2x5();
      photonB->e3x3               = it->e3x3();
      photonB->e5x5               = it->e5x5();
      photonB->r1x5               = it->r1x5();
      photonB->r2x5               = it->r2x5();
      photonB->maxEnergyXtal      = it->maxEnergyXtal();

      photonB->hasConversionTracks = it->hasConversionTracks();      
      if (it->hasConversionTracks()) {
        const reco::ConversionRefVector conversions = it->conversions();
        for (edm::RefVector<reco::ConversionCollection>::const_iterator jt  = conversions.begin();
                                                                        jt != conversions.end(); 
                                                                      ++jt) 
	{
          const reco::Conversion& obj = (**jt);
	  if (obj.nTracks() < 2 or
              !obj.conversionVertex().isValid()) continue;
          photonB->nTracks = obj.nTracks();
          photonB->isConverted = obj.isConverted();
          photonB->pairInvMass = obj.pairInvariantMass();
          photonB->pairCotThetaSeparation
     	        = obj.pairCotThetaSeparation();

	  math::XYZVectorF  mom = obj.pairMomentum();
          photonB->pairPx = mom.x();
          photonB->pairPy = mom.y();
          photonB->pairPz = mom.z();

          const reco::Vertex &vertex = obj.conversionVertex();
          photonB->conv_vx = vertex.x();
          photonB->conv_vy = vertex.y();
          photonB->conv_vz = vertex.z();

	  photonB->eovp              = obj.EoverP();
	  photonB->zpv               = obj.zOfPrimaryVertexFromTracks();
	  photonB->distOfMinApproach = obj.distOfMinimumApproach();
	  photonB->dPhiTracksAtVtx   = obj.dPhiTracksAtVtx();
	  photonB->dPhiTracksAtEcal  = obj.dPhiTracksAtEcal();
	  photonB->dEtaTracksAtEcal  = obj.dEtaTracksAtEcal();
        }    
      }
    }
  }
  else {
    edm::LogError("PhotonBlock") << "Error >> Failed to get pat::Photon for label: " 
                                 << _photonInputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhotonBlock);
