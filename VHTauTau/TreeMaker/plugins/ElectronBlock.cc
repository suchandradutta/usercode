#include <iostream>
#include <algorithm>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TClonesArray.h"

#include "VHTauTau/TreeMaker/plugins/ElectronBlock.h"
#include "VHTauTau/TreeMaker/interface/PhysicsObjects.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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
ElectronBlock::ElectronBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _trkInputTag(iConfig.getParameter<edm::InputTag>("trackSrc")),
  _dcsInputTag(iConfig.getParameter<edm::InputTag>("dcsSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  _electronInputTag(iConfig.getParameter<edm::InputTag>("electronSrc")),
  _pfElectronInputTag(iConfig.getParameter<edm::InputTag>("pfElectronSrc")),
  _ecalEBInputTag(iConfig.getParameter<edm::InputTag>("ecalEBInputTag")),
  _ecalEEInputTag(iConfig.getParameter<edm::InputTag>("ecalEEInputTag"))
{}
void ElectronBlock::beginJob() 
{
  // Get TTree pointer
  //  edm::Service<TFileService> fs;
  //fs->file().cd("/");
  //fs->cd();
  //  TFile& file = fs->file();
  //  TTree* tree = dynamic_cast<TTree*>(file.Get("vhtree")); 
  //assert(tree);

  std::string tree_name = "vhtree";
  TTree* tree = Utility::getTree(tree_name);
  cloneElectron = new TClonesArray("Electron");
  tree->Branch("Electron", &cloneElectron, 32000, 2);
  tree->Branch("nElectron", &fnElectron,  "fnElectron/I");
}
void ElectronBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneElectron->Clear();
  fnElectron = 0;

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(_trkInputTag, tracks);

  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByLabel(_dcsInputTag, dcsHandle);

#if 0
  edm::ESHandle<CaloGeometry> caloGeometry;
  iSetup.get<CaloGeometryRecord>().get(caloGeometry);

  edm::ESHandle<CaloTopology> caloTopology;
  iSetup.get<CaloTopologyRecord>().get(caloTopology);

  edm::Handle<EcalRecHitCollection> ecalBarrelRecHitHandle;
  iEvent.getByLabel(_ecalEBInputTag, ecalBarrelRecHitHandle);

  edm::Handle<EcalRecHitCollection> ecalEndcapRecHitHandle;
  iEvent.getByLabel(_ecalEEInputTag, ecalEndcapRecHitHandle);

  //  edm::Handle<reco::TrackCollection> trackHandle;
  //iEvent.getByLabel(trkInputTag, trackHandle);
  //const reco::TrackCollection* trackColl = tracks.product();
  edm::Handle<reco::BeamSpot> BeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot",BeamSpotHandle);
  const reco::BeamSpot* spot = BeamSpotHandle.product();
  PhotonTkIsolation TrackTool(0.3, 0.04, 0.7, 0.2, 9999, tracks.product(), math::XYZPoint(spot->x0(),spot->y0(),spot->z0()));

  EcalClusterLazyTools EcalTool(iEvent,iSetup,_ecalEBInputTag,_ecalEEInputTag);

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
  EcalRecHitMetaCollection ecalBarrelHits(*ecalBarrelRecHitHandle); // these are all needed to make an instance of EgammaRecHitIsolation
  EgammaRecHitIsolation ecalBarrelIsol(ConeOutRadius, ConeInRadius, EtaWidth, PtMin, EMin, 
                                            caloGeometry, &ecalBarrelHits, sevLevel, DetId::Ecal);
  EgammaRecHitIsolation HeepEcalBarrelIsol(HeepConeOutRadius, HeepConeInRadius, HeepEtaWidth, 
                                            PtMin, EMin, caloGeometry, &ecalBarrelHits, sevLevel, DetId::Ecal);
  HeepEcalBarrelIsol.setUseNumCrystals(true);
  EcalRecHitMetaCollection ecalEndcapHits(*ecalEndcapRecHitHandle);
  EgammaRecHitIsolation ecalEndcapIsol(ConeOutRadius, EndcapConeInRadius, EtaWidth, EndcapPtMin, EndcapEMin, 
                                            caloGeometry, &ecalEndcapHits, sevLevel, DetId::Ecal);
  EgammaRecHitIsolation HeepEcalEndcapIsol(HeepConeOutRadius, HeepConeInRadius, HeepEtaWidth, EndcapPtMin, 
                                            EndcapEMin, caloGeometry, &ecalEndcapHits, sevLevel, DetId::Ecal);
  HeepEcalEndcapIsol.setUseNumCrystals(true);
#endif
  double evt_bField = 3.8;
  // need the magnetic field
  //
  // if isRealData then derive bfield using the
  // magnet current from DcsStatus
  // otherwise take it from the IdealMagneticFieldRecord
  if (iEvent.isRealData()) {
    if (dcsHandle.isValid()) {
      edm::LogInfo("ElectronBlock") << "Successfully obtained " << _dcsInputTag;
      // scale factor = 3.801/18166.0 which are
      // average values taken over a stable two-week period
      double currentToBFieldScaleFactor = 2.09237036221512717e-04;
      double current = (*dcsHandle)[0].magnetCurrent();
      evt_bField = current*currentToBFieldScaleFactor;
    } 
    else {
      edm::LogError("ElectronBlock") << "Error! Can't get the product " << _dcsInputTag;
    }
  } 
  else {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
    if (magneticField.isValid()) {
      edm::LogInfo("ElectronBlock") << "Successfully obtained IdealMagneticFieldRecord!";
      evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
    } 
    else {
      edm::LogError("ElectronBlock") << "Error! Can't get IdealMagneticFieldRecord";
    }
  }

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(_electronInputTag, electrons);

  edm::Handle<reco::PFCandidateCollection> pfElectrons;
  iEvent.getByLabel(_pfElectronInputTag, pfElectrons);

  if (electrons.isValid()) {
    edm::LogInfo("ElectronBlock") << "Total # PAT Electrons: " << electrons->size();
    for (std::vector<pat::Electron>::const_iterator it = electrons->begin(); 
                                                   it != electrons->end(); ++it) {
      if (fnElectron == kMaxElectron) {
	edm::LogInfo("ElectronBlock") << "Too many PAT Electrons, fnElectron = " << fnElectron; 
	break;
      }
      // if electron is not ECAL driven, continue
      if (!it->ecalDrivenSeed()) continue;

      // Conversion
      ConversionFinder convFinder;
      double dist = -9999.;
      double dcot = -9999.;
      if (tracks.isValid()) {
	edm::LogInfo("ElectronBlock") << "Successfully obtained " << _trkInputTag;

        ConversionInfo convInfo = convFinder.getConversionInfo(*it, tracks, evt_bField);
        dist = convInfo.dist();
        dcot = convInfo.dcot();
      } 
      else {
	edm::LogError("ElectronBlock") << "Error! Can't get the product " << _trkInputTag;
      }

      // Vertex association
      double minVtxDist3D = 9999.;
      int indexVtx = -1;
      double vertexDistZ = 9999.;
 
      if (primaryVertices.isValid()) {
	edm::LogInfo("ElectronBlock") << "Total # Primary Vertices: " << primaryVertices->size();
        for (reco::VertexCollection::const_iterator v_it = primaryVertices->begin() ; 
                                                   v_it != primaryVertices->end() ; ++v_it) {
          double dist3D = sqrt(pow(it->gsfTrack()->dxy(v_it->position()),2) + pow(it->gsfTrack()->dz(v_it->position()),2));
          if (dist3D<minVtxDist3D) {
            minVtxDist3D = dist3D;
            indexVtx = int(std::distance(primaryVertices->begin(),v_it));
            vertexDistZ = it->gsfTrack()->dz(v_it->position());
          }
        }
      } 
      else {
	edm::LogError("ElectronBlock") << "Error! Can't get the product " << _vtxInputTag;
      }      
      double pfreliso = it->userFloat("pfLooseIsoPt04")/it->pt();

      electronB = new ((*cloneElectron)[fnElectron++]) Electron();
      electronB->eta        = it->eta();
      electronB->phi        = it->phi();
      electronB->pt         = it->pt();
      electronB->trackPt    = it->gsfTrack()->pt();
      electronB->energy     = it->energy();
      electronB->caloEnergy = it->caloEnergy();
      electronB->charge     = it->charge();

      // ID variables
      electronB->hoe           = it->hadronicOverEm();
      electronB->sigmaEtaEta   = it->sigmaEtaEta();
      electronB->sigmaIEtaIEta = it->sigmaIetaIeta();
      electronB->deltaPhiTrkSC = it->deltaPhiSuperClusterTrackAtVtx();
      electronB->deltaEtaTrkSC = it->deltaEtaSuperClusterTrackAtVtx();
      electronB->classif       = it->classification();
      electronB->e1x5overe5x5  = (it->e5x5() > 0) ? (it->e1x5()/it->e5x5()) : 0;
      electronB->e2x5overe5x5  = (it->e5x5() > 0) ? (it->e2x5Max()/it->e5x5()) : 0;

      // Iso variables
      electronB->isoEcal03 = it->dr03EcalRecHitSumEt();
      electronB->isoHcal03 = it->dr03HcalTowerSumEt();
      electronB->isoTrk03  = it->dr03TkSumPt();
      electronB->isoEcal04 = it->dr04EcalRecHitSumEt();
      electronB->isoHcal04 = it->dr04HcalTowerSumEt();
      electronB->isoTrk04  = it->dr04TkSumPt();
      electronB->isoRel03  = (it->dr03EcalRecHitSumEt()+it->dr03HcalTowerSumEt()+it->dr03TkSumPt())/it->pt();
      electronB->isoRel04  = (it->dr04EcalRecHitSumEt()+it->dr04HcalTowerSumEt()+it->dr04TkSumPt())/it->pt();

      // Conversion variables
      electronB->missingHits = it->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      electronB->dist_vec    = dist;
      electronB->dCotTheta   = dcot;

      // SC associated with electron
      electronB->scEn        = it->superCluster()->energy();
      electronB->scEta       = it->superCluster()->eta();
      electronB->scPhi       = it->superCluster()->phi();
      electronB->scET        = it->superCluster()->energy()/cosh(it->superCluster()->eta());
      electronB->scRawEnergy = it->superCluster()->rawEnergy();

      // Vertex association variables
      electronB->vtxDist3D = minVtxDist3D;
      electronB->vtxIndex  = indexVtx;
      electronB->vtxDistZ  = vertexDistZ;
      electronB->pfRelIso  = pfreliso;
#if 0
      // Ecal Spike Cleaning
      double emax    = -1.;
      double e9      = 1/999.;
      double eright  = 999.;
      double eleft   = 0.;
      double etop    = 0.;
      double ebottom = 0.;
      int ekoutoftime = -1;

      if (it->superCluster()->seed()->seed().subdetId() == EcalBarrel) {
        emax    = EcalClusterTools::eMax( *(it->superCluster()), &(*ecalBarrelRecHitHandle) );
        e9      = EcalClusterTools::e3x3( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) );
        eright  = EcalClusterTools::eRight( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) );
        eleft   = EcalClusterTools::eLeft( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) ) ;
        etop    = EcalClusterTools::eTop( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) ) ;
        ebottom = EcalClusterTools::eBottom( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) );

        int flaggedRecHitCounter = 0;
        const std::vector<std::pair<DetId, float> > & hitsAndFractions = it->superCluster()->seed()->hitsAndFractions();
	std::vector<std::pair<DetId, float> >::const_iterator hitIter;
        for (hitIter = hitsAndFractions.begin(); hitIter != hitsAndFractions.end(); ++hitIter) {
	  EcalRecHitCollection::const_iterator recHit = ecalBarrelRecHitHandle->find(hitIter->first);
          if (recHit != ecalBarrelRecHitHandle->end()) {
            if ( (hitIter->second*recHit->energy()/it->superCluster()->rawEnergy() > 0.05)
		&& recHit->checkFlag(EcalRecHit::kOutOfTime) ) flaggedRecHitCounter++;
          }
        }
        ekoutoftime = (flaggedRecHitCounter > 0) ? 1 : 0;
      } 
      else {
        emax    = EcalClusterTools::eMax( *(it->superCluster()), &(*ecalEndcapRecHitHandle) );
        e9      = EcalClusterTools::e3x3( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) );
        eright  = EcalClusterTools::eRight( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) );
        eleft   = EcalClusterTools::eLeft( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) ) ;
        etop    = EcalClusterTools::eTop( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) ) ;
        ebottom = EcalClusterTools::eBottom( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) );

        int flaggedRecHitCounter = 0;
        const std::vector<std::pair<DetId, float> > & hitsAndFractions = it->superCluster()->seed()->hitsAndFractions();
	std::vector<std::pair<DetId, float> >::const_iterator hitIter;
        for (hitIter = hitsAndFractions.begin(); hitIter != hitsAndFractions.end(); ++hitIter) {
	  EcalRecHitCollection::const_iterator recHit = ecalEndcapRecHitHandle->find(hitIter->first);
          if (recHit != ecalEndcapRecHitHandle->end()) {
            if ( (hitIter->second*recHit->energy()/it->superCluster()->rawEnergy()) > 0.05 
		 && recHit->checkFlag(EcalRecHit::kOutOfTime) ) flaggedRecHitCounter++;
          }
        }
        ekoutoftime = (flaggedRecHitCounter>0) ? 1 : 0;
      }

      electronB->scE1E9 = emax/e9;
      electronB->scS4S1 = (eright+eleft+etop+ebottom)/emax; // 
      electronB->sckOutOfTime = ekoutoftime;

      reco::SuperClusterRef eleSCRef = it->superCluster();  // get SCRef to use to make ele candidate
      TVector3 sc_vec;
      sc_vec.SetXYZ(eleSCRef->x(),eleSCRef->y(),eleSCRef->z());
      double eleSCRefpt = eleSCRef->energy()*(sc_vec.Perp()/sc_vec.Mag());
      reco::RecoEcalCandidate ecalCand;  // make ele candidate to use Iso algorithm
      ecalCand.setSuperCluster(eleSCRef);
      const reco::Candidate::PolarLorentzVector photon_vec(eleSCRefpt, eleSCRef->eta(), eleSCRef->phi(), 0.0);
      ecalCand.setP4(photon_vec);

      electronB->scEcalIso = (fabs(eleSCRef->eta()) < 1.48) ? ecalBarrelIsol.getEtSum(&ecalCand)
                                                            : ecalEndcapIsol.getEtSum(&ecalCand);
      electronB->scHEEPEcalIso = (fabs(eleSCRef->eta()) < 1.48) ? HeepEcalBarrelIsol.getEtSum(&ecalCand)
                                                                : HeepEcalEndcapIsol.getEtSum(&ecalCand);
      electronB->scHEEPTrkIso = TrackTool.getPtTracks(&ecalCand);
#endif
    }
  } 
  else {
    edm::LogError("ElectronBlock") << "Error! Can't get the product " << _electronInputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronBlock);
