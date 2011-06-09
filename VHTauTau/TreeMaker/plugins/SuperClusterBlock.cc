#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/SuperClusterBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/PhotonTkIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

#include "TVector3.h"

// Constructor
SuperClusterBlock::SuperClusterBlock(const edm::ParameterSet& iConfig) :
  _tree(0), 
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _ebInputTag(iConfig.getParameter<edm::InputTag>("EBInputTag")),
  _eeInputTag(iConfig.getParameter<edm::InputTag>("EEInputTag")),
  _ecalEBInputTag(iConfig.getParameter<edm::InputTag>("EcalEBInputTag")),
  _ecalEEInputTag(iConfig.getParameter<edm::InputTag>("EcalEEInputTag")),
  _trkInputTag(iConfig.getParameter<edm::InputTag>("TracksInputTag")),
  _eleInputTag(iConfig.getParameter<edm::InputTag>("ElectronsInputTag"))
{}
void SuperClusterBlock::beginJob() 
{
  std::string tree_name = "vhtree";
  if (_tree) _tree = Utility::getTree(tree_name);
  cloneSuperCluster = new TClonesArray("SuperCluster");
  _tree->Branch("SuperCluster", &cloneSuperCluster, 32000, 2);
  _tree->Branch("nSuperCluster", &fnSuperCluster,  "fnSuperCluster/I");
}
void SuperClusterBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneSuperCluster->Clear();
  fnSuperCluster = 0;

  edm::Handle<reco::SuperClusterCollection> superClustersEBHandle;
  iEvent.getByLabel(_ebInputTag, superClustersEBHandle);

  edm::Handle<reco::SuperClusterCollection> superClustersEEHandle;
  iEvent.getByLabel(_eeInputTag, superClustersEEHandle);

  edm::ESHandle<CaloGeometry> caloGeometry;
  iSetup.get<CaloGeometryRecord>().get(caloGeometry);

  edm::ESHandle<CaloTopology> caloTopology;
  iSetup.get<CaloTopologyRecord>().get(caloTopology);

  edm::Handle<EcalRecHitCollection> ecalBarrelRecHitHandle;
  iEvent.getByLabel(_ecalEBInputTag, ecalBarrelRecHitHandle);

  edm::Handle<EcalRecHitCollection> ecalEndcapRecHitHandle;
  iEvent.getByLabel(_ecalEEInputTag, ecalEndcapRecHitHandle);

  EcalClusterLazyTools EcalTool(iEvent, iSetup, _ecalEBInputTag, _ecalEEInputTag);

  edm::Handle<reco::TrackCollection> trackHandle;
  iEvent.getByLabel(_trkInputTag, trackHandle);
  const reco::TrackCollection* tracks = trackHandle.product();
  edm::Handle<reco::BeamSpot> BeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", BeamSpotHandle);
  const reco::BeamSpot* spot = BeamSpotHandle.product();
  PhotonTkIsolation TrackTool(0.3,0.04,0.7,0.2,9999,tracks,math::XYZPoint(spot->x0(),spot->y0(),spot->z0()));

  edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
  iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
  const EcalSeverityLevelAlgo* sevLevel = sevlv.product();

  double ConeOutRadius = 0.4;  // these are all needed to make an instance of EgammaRecHitIsolation
  double ConeInRadius = 0.045;
  double EtaWidth = 0.02;
  double PtMin = 0.;
  double EMin = 0.08;
  double EndcapConeInRadius = 0.07; // note that this is in number-of-crystal units
  double EndcapPtMin = 0.;
  double EndcapEMin = 0.1;
  double HeepConeOutRadius = 0.3;  // HEEP uses different iso values than Egamma/PAT default
  double HeepConeInRadius = 3.; // note that this is in num of crystals
  double HeepEtaWidth = 1.5;  // note that this is in num of crystals
  EcalRecHitMetaCollection ecalBarrelHits(*ecalBarrelRecHitHandle);// these are all needed to make an instance of EgammaRecHitIsolation
  EgammaRecHitIsolation ecalBarrelIsol(ConeOutRadius,ConeInRadius,EtaWidth,PtMin,EMin,caloGeometry,&ecalBarrelHits,sevLevel, DetId::Ecal);
  EgammaRecHitIsolation HeepEcalBarrelIsol(HeepConeOutRadius,HeepConeInRadius,HeepEtaWidth,PtMin,EMin,caloGeometry,&ecalBarrelHits,sevLevel, DetId::Ecal);
  HeepEcalBarrelIsol.setUseNumCrystals(true);
  EcalRecHitMetaCollection ecalEndcapHits(*ecalEndcapRecHitHandle);
  EgammaRecHitIsolation ecalEndcapIsol(ConeOutRadius,EndcapConeInRadius,EtaWidth,EndcapPtMin,EndcapEMin,caloGeometry,&ecalEndcapHits,sevLevel, DetId::Ecal);
  EgammaRecHitIsolation HeepEcalEndcapIsol(HeepConeOutRadius,HeepConeInRadius,HeepEtaWidth,EndcapPtMin,EndcapEMin,caloGeometry,&ecalEndcapHits,sevLevel, DetId::Ecal);
  HeepEcalEndcapIsol.setUseNumCrystals(true);

  // SuperClusters for barrel and endcap are in different collections
  // "hybrid" = barrel, "multi5x5" = endcap
  // Loop over both collections
  for (reco::SuperClusterCollection::const_iterator it = superClustersEBHandle->begin(); 
                                                   it != superClustersEBHandle->end(); ++it) {
    if (fnSuperCluster == kMaxSuperCluster) {
      edm::LogInfo("SuperClusterBlock") << "Too many EB SuperCluster, fnSuperCluster = " << fnSuperCluster; 
      break;
    }

    scB = new ((*cloneSuperCluster)[fnSuperCluster++]) SuperCluster();

    scB->eta = it->eta();
    scB->phi = it->phi();
    TVector3 sc_vec;
    sc_vec.SetXYZ(it->x(),it->y(),it->z());
    scB->pt = it->energy()*(sc_vec.Perp()/sc_vec.Mag());
    scB->rawEnergy = it->rawEnergy();
    scB->clustersSize = it->clustersSize();

    const reco::SuperCluster* pnt_sc = &(*it);
    HoECalculator calc_HoE; // requires HCAL RecHits that are only available in RECO files
    double schoe = calc_HoE(pnt_sc,iEvent,iSetup);
    scB->hoe = schoe;

    std::vector<float> scLocalCov = EcalTool.scLocalCovariances(*it);
    double scSigmaiEiE = sqrt(scLocalCov[0]); // same method used in GsfElectronAlgo.cc
    scB->sigmaIEtaIEta = scSigmaiEiE;

    reco::SuperClusterRef tempSCRef(superClustersEBHandle,std::distance(superClustersEBHandle->begin(),it));  // get SCRef to use to make ele candidate
    reco::RecoEcalCandidate ecalCand;  // make ele candidate to use Iso algorithm
    ecalCand.setSuperCluster(tempSCRef);
    scB->ecalIso = ecalBarrelIsol.getEtSum(&ecalCand);
    scB->heepTrkIso = TrackTool.getPtTracks(&ecalCand);
    scB->heepEcalIso = HeepEcalBarrelIsol.getEtSum(&ecalCand);

    double closestTrkDr = 99.;
    reco::Track closestTrk;
    double nextTrkDr = 99.;
    reco::Track nextTrk;
    for (reco::TrackCollection::const_iterator trk = trackHandle->begin(); 
                                              trk != trackHandle->end(); ++trk) {
      TVector3 trk_vec;
      trk_vec.SetPtEtaPhi(trk->pt(), trk->eta(), trk->phi());
      double dR = trk_vec.DeltaR(sc_vec);
      if (dR < closestTrkDr) {
        nextTrkDr = closestTrkDr;
        nextTrk = closestTrk;
        closestTrkDr = dR;
        closestTrk = *trk;
      } 
      else if (dR < nextTrkDr) {
        nextTrkDr = dR;
        nextTrk = *trk;
      }
    }
    scB->trackMatch = (closestTrkDr < 0.04) ? 1 : 0;  // 1=true, 0=false
    scB->dRTrack1 = closestTrkDr;
    scB->dRTrack2 = nextTrkDr;
    scB->track1Eta = closestTrk.eta();
    scB->track2Eta = nextTrk.eta();
    scB->track1Phi = closestTrk.phi();
    scB->track2Phi = nextTrk.phi();
    scB->track1Pt = closestTrk.pt();
    scB->track2Pt = nextTrk.pt();

    double emax    = EcalClusterTools::eMax( *it, &(*ecalBarrelRecHitHandle) );
    double e9      = EcalClusterTools::e3x3( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) );
    double eright  = EcalClusterTools::eRight( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) );
    double eleft   = EcalClusterTools::eLeft( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) ) ;
    double etop    = EcalClusterTools::eTop( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) ) ;
    double ebottom = EcalClusterTools::eBottom( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) );

    scB->e1e9 = emax/e9;
    scB->s4s1 = (eright+eleft+etop+ebottom)/emax;

    int flaggedRecHitCounter = 0;
    const std::vector<std::pair<DetId, float> > & hitsAndFractions = it->seed()->hitsAndFractions();
    std::vector<std::pair<DetId, float> >::const_iterator hitIter;

    for (hitIter = hitsAndFractions.begin(); hitIter != hitsAndFractions.end(); ++hitIter) {
      EcalRecHitCollection::const_iterator recHit = ecalBarrelRecHitHandle->find(hitIter->first);
      if (recHit != ecalEndcapRecHitHandle->end()) {
        if (hitIter->second*recHit->energy()/it->rawEnergy() > 0.05 
          && recHit->checkFlag(EcalRecHit::kOutOfTime)) flaggedRecHitCounter++;
      }
    }

    scB->kOutOfTime = (flaggedRecHitCounter > 0) ? 1 : 0;
  }
  for (reco::SuperClusterCollection::const_iterator it = superClustersEEHandle->begin(); 
                                                   it != superClustersEEHandle->end(); ++it) {
    if (fnSuperCluster == kMaxSuperCluster) {
      edm::LogInfo("SuperClusterBlock") << "Too many EE SuperCluster, fnSuperCluster = " << fnSuperCluster; 
      break;
    }
    scB = new ((*cloneSuperCluster)[fnSuperCluster++]) SuperCluster();
    scB->eta = it->eta();
    scB->phi = it->phi();
    TVector3 sc_vec;
    sc_vec.SetXYZ(it->x(), it->y(), it->z());
    scB->pt = it->energy()*(sc_vec.Perp()/sc_vec.Mag());
    scB->rawEnergy = it->rawEnergy();
    scB->clustersSize = it->clustersSize();

    const reco::SuperCluster* pnt_sc = &(*it);
    HoECalculator calc_HoE; // requires HCAL RecHits that are only available in RECO files
    double schoe = calc_HoE(pnt_sc,iEvent,iSetup);
    scB->hoe = schoe;

    std::vector<float> scLocalCov = EcalTool.scLocalCovariances(*it);
    double scSigmaiEiE = sqrt(scLocalCov[0]); // same method used in GsfElectronAlgo.cc
    scB->sigmaIEtaIEta = scSigmaiEiE;

    reco::SuperClusterRef tempSCRef(superClustersEEHandle,std::distance(superClustersEEHandle->begin(),it));  // get SCRef to use to make ele candidate
    reco::RecoEcalCandidate ecalCand;  // make ele candidate to use Iso algorithm
    ecalCand.setSuperCluster(tempSCRef);
    scB->ecalIso = ecalEndcapIsol.getEtSum(&ecalCand);
    scB->heepEcalIso = HeepEcalEndcapIsol.getEtSum(&ecalCand);
    scB->heepTrkIso = TrackTool.getPtTracks(&ecalCand);

    double closestTrkDr = 99.;
    reco::Track closestTrk;
    double nextTrkDr = 99.;
    reco::Track nextTrk;
    for (reco::TrackCollection::const_iterator trk = trackHandle->begin(); 
                                              trk != trackHandle->end(); ++trk) {
      TVector3 trk_vec;
      trk_vec.SetPtEtaPhi(trk->pt(), trk->eta(), trk->phi());
      double dR = trk_vec.DeltaR(sc_vec);
      if (dR < closestTrkDr) {
        nextTrkDr = closestTrkDr;
        nextTrk = closestTrk;
        closestTrkDr = dR;
        closestTrk = *trk;
      } 
      else if (dR < nextTrkDr) {
        nextTrkDr = dR;
        nextTrk = *trk;
      }
    }
    scB->trackMatch = (closestTrkDr < 0.04) ? 1 : 0;  // 1=true, 0=false
    scB->dRTrack1 = closestTrkDr;
    scB->dRTrack2 = nextTrkDr;
    scB->track1Eta = closestTrk.eta();
    scB->track2Eta = nextTrk.eta();
    scB->track1Phi = closestTrk.phi();
    scB->track2Phi = nextTrk.phi();
    scB->track1Pt = closestTrk.pt();
    scB->track2Pt = nextTrk.pt();

    double emax    = EcalClusterTools::eMax( *it, &(*ecalEndcapRecHitHandle) );
    double e9      = EcalClusterTools::e3x3( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) );
    double eright  = EcalClusterTools::eRight( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) );
    double eleft   = EcalClusterTools::eLeft( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) ) ;
    double etop    = EcalClusterTools::eTop( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) ) ;
    double ebottom = EcalClusterTools::eBottom( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) );

    scB->e1e9 = emax/e9;
    scB->s4s1 = (eright+eleft+etop+ebottom)/emax;

    int flaggedRecHitCounter = 0;
    const std::vector<std::pair<DetId, float> > & hitsAndFractions = it->seed()->hitsAndFractions();
    std::vector<std::pair<DetId, float> >::const_iterator hitIter;
    for (hitIter=hitsAndFractions.begin(); hitIter!=hitsAndFractions.end(); ++hitIter ) {
      EcalRecHitCollection::const_iterator recHit = ecalEndcapRecHitHandle->find(hitIter->first);
      if (recHit != ecalEndcapRecHitHandle->end()) {
        if (hitIter->second*recHit->energy()/it->rawEnergy() > 0.05 
	     && recHit->checkFlag(EcalRecHit::kOutOfTime)) flaggedRecHitCounter++;
      }
    }
    scB->kOutOfTime = (flaggedRecHitCounter>0) ? 1 : 0;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SuperClusterBlock);
