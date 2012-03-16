#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "MuTauTau.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::abs;
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
MuTauTau::MuTauTau()
  : AnaBase()
{}
// ----------
// Destructor
// ----------
MuTauTau::~MuTauTau() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MuTauTau::beginJob() 
{ 
  AnaBase::beginJob();

  count_sig_event = 0;
  for (unsigned int i = 0; i < NEL(nEvtSel); ++i) {
    nEvtSel[i] = 0;
  }
  _histf->cd();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void MuTauTau::bookHistograms() 
{
  static const double pi = TMath::Pi();
  new TH1F("nMuon", "Number of Muons", 6, -0.5, 5.5);
  new TH1F("muonPt1_allcut", "Highest Pt Muon Pt Distribution after all cut applied", 100, -0.5, 99.5);
  new TH1F("muonPt1_stage1", "Highest Pt Muon Pt Distribution after vtx+Muon eta ID ISO and Pt cut", 100, -0.5, 99.5);
  new TH1F("muonPt1_stage2", "Highest Pt Muon Pt Distribution after vtx+Muon eta ID ISO and Pt and 1 Tau with Pt cut", 100, -0.5, 99.5);
  new TH1F("muonPt1_stage3", "Highest Pt Muon Pt Distribution after vtx+Muon eta ID ISO and Pt and 2 Tau with Pt cut", 100, -0.5, 99.5);
  new TH1F("muonEta1_allcut", "Highest Pt Muon Eta Distribution after all cut applied", 100, -3, 3);
  new TH1F("muonEta1_stage1", "Highest Pt Muon Eta Distribution after vtx+Muon eta ID ISO and Pt cut", 100, -3, 3);
  new TH1F("muonEta1_stage2", "Highest Pt Muon Eta Distribution after vtx+Muon eta ID ISO and Pt and 1 Tau with Pt cut", 100, -3, 3);
  new TH1F("muonEta1_stage3", "Highest Pt Muon Eta Distribution after vtx+Muon eta ID ISO and Pt and 2 Tau with Pt cut", 100, -3, 3);
  new TH1F("tauPt1_allcut", "Highest Pt Tau Pt Distribution after all cut applied", 100, -0.5, 99.5);
  new TH1F("tauPt1_stage1", "Highest Pt Tau Pt Distribution after vtx+Muon eta ID ISO and Pt and 1 Tau with Pt cut", 100, -0.5, 99.5);
  new TH1F("tauPt1_stage2", "Highest Pt Tau Pt Distribution after vtx+Muon eta ID ISO and Pt and 2 Tau with Pt cut", 100, -0.5, 99.5);
  new TH1F("tauEta1_allcut", "Highest Pt Tau Eta Distribution after all cut applied", 100, -3, 3);
  new TH1F("tauEta1_stage1", "Highest Pt Tau Eta Distribution after vtx+Muon eta ID ISO and Pt and 1 Tau with Pt cut", 100, -3, 3);
  new TH1F("tauEta1_stage2", "Highest Pt Tau Eta Distribution after vtx+Muon eta ID ISO and Pt and 2 Tau with Pt cut", 100, -3, 3);
  new TH1F("tauPt2_allcut", "2nd Highest Pt Tau Pt Distribution after all cut applied", 100, -0.5, 99.5);
  new TH1F("tauPt2_stage1", "2nd Highest Pt Tau Pt Distribution after vtx+Muon eta ID ISO and Pt and 2 Tau with Pt cut", 100, -0.5, 99.5);
  new TH1F("tauEta2_allcut", "2nd Highest Pt Tau Eta Distribution after all cut applied", 100, -3, 3);
  new TH1F("tauEta2_stage1", "2nd Highest Pt Tau Eta Distribution after vtx+Muon eta ID ISO and Pt and 2 Tau with Pt cut", 100, -3, 3);
  new TH1F("mumuInv","Invariant mass of two muons if there exists more than 1 good muon", 140, 0, 140);
  new TH1F("tautauInvmass", "Tau Tau Invariant mass// vissible Higgs mass", 150, -0.5, 149.5);
  new TH1F("muonPt2", "Second highest Pt Muon Pt Distribution", 100, -0.5, 99.5);
  new TH1F("nTau", "Number of Taus", 10, -0.5, 9.5);
  new TH1F("nTau_gt5", "Number of Taus with pt more than 5", 10, -0.5, 9.5);
  new TH1F("tauPt3", "3rd Tau Pt Distribution", 100, -0.5, 99.5);
  new TH1F("deltaR_mutau1", "delta R difference between the selected muon and highest pt tau after all selection requirment MUTAUTAU channel", 100, 0, 10);
  new TH1F("deltaR_mutau2", "delta R difference between the selected muon and 2nd highest pt tau after all selection requirment MUTAUTAU channel", 100, 0, 10);
  new TH1F("deltaR_tau1tau2", "delta R difference between the selected two tau's after all selection requirment MUTAUTAU channel", 100, 0, 10);
  new TH1F("met", "met distribution of the event", 140, 0, 140); 
  new TH1F("MT_Mu_Met","transverse mass plot for Muon and MET", 140, 0, 140);
  new TH1F("MT_Tau1_Met","transverse mass plot for Tau1 and MET", 140, 0, 140);
  new TH1F("MT_Tau2_Met","transverse mass plot for Tau2 and MET", 140, 0, 140);
  new TH1F("DiTauPt", "Pt Distribution of the Di-Tau ", 140, 0, 140);
  new TH1F("var1", " plot of  Pt(Tau1+Tau2)/(tauiPt + tau2Pt)", 100, 0, 2);
  new TH1F("diffPhi_Mu_Met", "delta Phi difference of Mu and Met", 100, 0, 5);
  new TH1F("diffPhi_Mu_DiTau", "delta Phi difference of Mu and Di-Tau", 100, 0, 5);
  new TH1F("diffPhi_Mu_Tau1", "delta Phi difference of Mu and Tau1", 100, 0, 5);

  //new TH1F("deltaPt1","minimum Pt difference  between detector and generator level muon1", 100, 0, 100);
  //new TH1F("deltaPt2","minimum Pt difference  between detector and generator level muon2", 100, 0, 100);
  //new TH1F("deltaEta1","minimum Eta difference  between detector and generator level muon1", 100, 0, 100);
  //new TH1F("deltaEta2","minimum Eta difference  between detector and generator level muon2", 100, 0, 100);
  //new TH1F("deltaPhi1","minimum Phi difference  between detector and generator level muon1", 100, 0, 100);
  //new TH1F("deltaPhi2","minimum Phi difference  between detector and generator level muon2", 100, 0, 100);
  //new TH1F("deltaDR","minimum DR difference  between detector and generator level muon", 100, 0, 100);
  new TH1F("genInvMass", "Invariant mass of two muons in generator level irrespective their mother", 100, 0, 140);
  new TH1F("oddInvMass", "Inv Mass of the two muons in generator level requiring that the mothers of two muons are either W or Tau", 100, 0, 140);
  
  new TH1F("PhiMu1","phi angle of heighst pt muon", 100, -7, 7);
  new TH1F("PhiMu2","phi angle of 2nd heighst pt muon", 100, -7, 7);
  new TH1F("PhiTau","phi angle of heighst pt Tau", 100, -7, 7);

  new TH1F("WH_DR", "DR difference between W and Higgs at the generator level", 100, 0, 10);
  new TH1F("Ht1_DR","DR diff between the Higgs and the Hadronic Tau at the generator level", 100, 0, 10);
  new TH1F("HM2_DR","DR diff betwen the Higgs and tau decayed Muon at the generator level", 100, 0, 10);
  new TH1F("WM1_DR","DR diff between the W and first muon at the generator level", 100, 0, 10);
  new TH1F("t1M2_DR","DR diff between the Hadronic Tau and Tau decayed Muon at the generator level", 100, 0, 10);
  new TProfile("WH_Wpt","DR diff between W and H as a function of W's pt", 100, 0, 100, 0, 10);
  new TProfile("t1M2_Hpt","DR diff between Hadronic Tau and Tau decayed Muon as a function of Higgs's pt",100, 0, 100, 0, 10);
  new TProfile("WM1_Wpt", "DR diff between W boson and daughter Muon as a function of W's pt", 100, 0, 100, 0,10);
  new TH1F("DPhi_WH", "Phi diff between W and H", 100, -10, 10);
  new TH1F("DEta_WH", "Eta diff between W and H", 100, -10, 10);
  new TH1F("DPhi_t1M2","Phi diff between Hadronic Tau and Tau decayed Muon", 100, -10, 10);
  new TH1F("DEta_t1M2","Eta diff between Hadronic Tau and Tau decayed Muon", 100, -10, 10);
  new TH1F("mumuInvMass_ZH", "inv mass of two muons coming from Z decay for ZH process", 100, 0, 140);

  if (_isMC) {
    new TH1F("npu", "nPileUp", 50, 0, 50);
    new TH1F("puwt", "PileUp weight factor", 100, 0, 2.0);
  } 
  new TH1I("ecounter", "Selected event counter", 15, -0.5, 14.5);
  new TH1F("nvertex", "No of good vertex", 30, -0.5, 29.5);
}
// -------------------
// The main event loop
// -------------------
void MuTauTau::clearLists() {
  vtxList.clear();
  muoList.clear();
  eleList.clear();
  tauList.clear();
  bjetList.clear();

  genMuonList.clear();
  genTauList.clear();
  GEN_MUON_LIST.clear();
  genMuoaList.clear();
  genMuobList.clear();
  genHList.clear();
  genWList.clear();
  genTau_H_List.clear();
  genMuon_Z_List.clear();
}
void MuTauTau::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;
  
  int nPrint = max(10000, nEvents/1000);

  Options op;
  op.verbose = false;
  op.usesbit = true;

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  double nEvt = 0, nEvtTrig = 0, nEvtVtx = 0;
  for (int ev = 0; ev < nEvents; ++ev) {
    clearEvent();
    int lflag = _chain->LoadTree(ev); 
    int nbytes = getEntry(lflag);    // returns total bytes read

    // PileUP weight
    _puevWt = 1; // for data
    if (_isMC) {
      int npu = 0;
      _puevWt = wtPileUp(npu);
      AnaUtil::fillHist1D("npu", npu, 1.0);
      AnaUtil::fillHist1D("puwt", _puevWt, 1.0);
    }
    nEvt += _puevWt;
    AnaUtil::fillHist1D("ecounter", 0, _puevWt);
    
    // Trigger selection
    if (_studyTrigger && !isTriggered()) continue;
    nEvtTrig += _puevWt;
    AnaUtil::fillHist1D("ecounter", 1, _puevWt);
    
    string currentFile(gSystem->BaseName(_chain->GetCurrentFile()->GetName())); 

    const Event* evt = dynamic_cast<Event*>(eventA->At(0));
    assert(evt);
    int run   = evt->run;
    int event = evt->event;

    // Show status of the run
    if (currentFile != lastFile) 
    cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
         << " ==> " << currentFile 
         << " <<< Run# " << run
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
    cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
         << " ==> " << _chain->GetCurrentFile()->GetName() 
         << " <<< Run# " << run 
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;

    if (_logOption > 0) 
    _fLog << "run: " << run
          << ", event: " << event
          << ", n_tau: "<< n_tau
          << ", n_muon: "<< n_muon
          << ", n_jet: " << n_jet
          << ", n_vertex: " << n_vertex
          << ", n_met: " << n_met
          << ", n_electron: " << n_electron 
          << endl;

    clearLists();
 
    //findGenInfo();
    if (_logOption >> 6 & 0x1) {
      bool dumpAll = false; // change to true to dump all the events
      if (dumpAll || (genMuonList.size() > 1 && genTauList.size() > 0)) dumpGenInfo(_fLog);
    }

    op.verbose = (_logOption >> 1 & 0x1); 
    findVtxInfo(vtxList, _vtxCutMap, op, _fLog);
    int nvtx = vtxList.size();
    double vz = ( nvtx > 0 ) ? vtxList.at(0).z : -999;
    AnaUtil::fillHist1D("nvertex", nvtx, _puevWt);

    op.verbose = (_logOption >> 2 & 0x1); 
    findTauInfo(tauList, _tauCutMap, vz, op, _fLog);

    op.verbose = (_logOption >> 3 & 0x1); 
    findMuonInfo(muoList, _muonCutMap, op, _fLog);

    op.verbose = (_logOption >> 4 & 0x1); 
    findElectronInfo(eleList, _electronCutMap, op, _fLog);

    op.verbose = (_logOption >> 5 & 0x1); 
    findJetInfo(bjetList, _bjetCutMap, op, _fLog);

    if (_logOption > 0)
    _fLog << "run: " << run
          << ", event: " << event
          << ", n_vertex_good: " << nvtx
          << ", n_muon_selected: " << muoList.size()
          << ", n_electron_selected: " << eleList.size()
          << ", n_tau_selected: " << tauList.size()
          << ", n_bjet_selected: " << bjetList.size()
          << endl;

    // Event Selection Starts here .....
    // presence of > 1 good vertex
    if (vtxList.size() < 1) continue;
    nEvtVtx += _puevWt;
    AnaUtil::fillHist1D("ecounter", 2, _puevWt);

    selectEvent();
    //computeDeltaR(muoList,tauList); 
  }  
  // Analysis is over
  nEvtSel[0] = nEvt;
  nEvtSel[1] = nEvtTrig;
  nEvtSel[2] = nEvtVtx;

  endJob();
}
void MuTauTau::findGenInfo() {
      // Generator level information
    if (n_genparticle) {
      //if (_logOption >> 6 & 0x1) 
      //  _fLog << "indx    status    pdgId     eta      phi      pt     energy            mID                             dID"
      //        << endl;
      for (int indx = 0; indx < n_genparticle; ++indx) {
        const GenParticle* gp = dynamic_cast<GenParticle*>(genParticleA->At(indx));
        if (!gp) continue;
        
	// cout << "indx = " << indx << ", id=" << gp->pdgId << ", status = " << gp->status << endl;
     
	//looking for a Hadronic Tau whose mother is Higgs
        if (gp->status == 2 && abs(gp->pdgId)==15) {
           vector<int> m = gp->motherIndices;
           const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(m[0]));
           if (!mgp) continue;

           int pdgid = mgp->pdgId;
           if (abs(pdgid) == 15) {
             vector<int> m2 = mgp->motherIndices;
             const GenParticle* m2gp = dynamic_cast<GenParticle*>(genParticleA->At(m2[0]));
             if (!m2gp) continue;
             pdgid = m2gp->pdgId;
             if (abs(pdgid)== 25) {
	       genHList.push_back(*m2gp);
             }
           }
           
           if (abs(pdgid) != 25) continue;

           vector<int> d = gp->daughtIndices;
           for (size_t i = 0; i < d.size(); ++i) {
	     int di = d[i];
             if (di >= n_genparticle) continue;
             const GenParticle* dgp = dynamic_cast<GenParticle*>(genParticleA->At(di));
             if (!dgp) continue;
             if (abs(dgp->pdgId) == 16) continue;
             if (abs(dgp->pdgId)!=11 && abs(dgp->pdgId)!=13 && abs(dgp->pdgId) !=12 && abs(dgp->pdgId) !=14 && dgp->pdgId !=22) {
               //cout << "pId=" << gp->pdgId << ", dId=" << dgp->pdgId << endl;
  	       genTauList.push_back(*gp);
               break;
             }
           }
          
        }


        // Looking for Muon whose mother is either W boson or Tau 
        if (gp->status == 1 && abs(gp->pdgId)== 13) {
	  vector<int> m = gp->motherIndices;
          const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(m[0]));
          if (!mgp) continue;
          if (abs (mgp->pdgId)==15)
	    genMuobList.push_back(*gp);
          if (abs (mgp->pdgId)==24){
	    genMuoaList.push_back(*gp);
            genWList.push_back(*mgp);
          }

          int mmid= mgp->pdgId;
          if (abs(mmid)==13){
            vector<int> m2 = mgp->motherIndices;
            const GenParticle* m2gp = dynamic_cast<GenParticle*>(genParticleA->At(m2[0]));
            if (!m2gp) continue;
            mmid = m2gp->pdgId;
            if (abs(mmid)==15)
	      genMuobList.push_back(*gp);
            if (abs(mmid)==24){
              genWList.push_back(*m2gp);
              genMuoaList.push_back(*gp);
	    }
          }
          
          
	   if (abs(mmid)==24 || abs(mmid)==15)
             genMuonList.push_back(*gp);
       	}

      	//looking for a Hadronic Tau whose mother is H-Boson
        if (gp->status == 2 && abs(gp->pdgId)==15) {
           vector<int> m = gp->motherIndices;
           const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(m[0]));
           if (!mgp) continue;

           int pdgid = mgp->pdgId;
           if (abs(pdgid) == 15) {
             vector<int> m2 = mgp->motherIndices;
             const GenParticle* m2gp = dynamic_cast<GenParticle*>(genParticleA->At(m2[0]));
             if (!m2gp) continue;
             pdgid = m2gp->pdgId;
           }
           
           if (abs(pdgid) != 25) continue;

           vector<int> d = gp->daughtIndices;
           for (size_t i = 0; i < d.size(); ++i) {
	     int di = d[i];
             if (di >= n_genparticle) continue;
             const GenParticle* dgp = dynamic_cast<GenParticle*>(genParticleA->At(di));
             if (!dgp) continue;
             if (abs(dgp->pdgId) == 16) continue;
             if (abs(dgp->pdgId)!=11 && abs(dgp->pdgId)!=13 && abs(dgp->pdgId) !=12 && abs(dgp->pdgId) !=14 && dgp->pdgId !=22) {
               //cout << "pId=" << gp->pdgId << ", dId=" << dgp->pdgId << endl;
  	       genTau_H_List.push_back(*gp);
               break;
             }
           }
          
        }
        // Looking for 2 Muon whose mother is Z 
        if (gp->status == 1 && abs(gp->pdgId)== 13) {
	  vector<int> m = gp->motherIndices;
          const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(m[0]));
          if (!mgp) continue;

          int mmid= mgp->pdgId;
          if (abs(mmid)==13){
            vector<int> m2 = mgp->motherIndices;
            const GenParticle* m2gp = dynamic_cast<GenParticle*>(genParticleA->At(m2[0]));
            if (!m2gp) continue;
            mmid = m2gp->pdgId;
            if (abs (mmid)== 13) {
              vector<int>m3 = m2gp->motherIndices;
              const GenParticle* m3gp = dynamic_cast<GenParticle*>(genParticleA->At(m3[0]));
              if (!m3gp) continue;
              mmid = m3gp->pdgId;
	    }
          }
 
	  if (abs(mmid)==23)
            genMuon_Z_List.push_back(*gp);
       	}



	//calculation of invariant mass in the generator level
	if(gp->status==1 && abs(gp->pdgId)== 13)
	  GEN_MUON_LIST.push_back(*gp); 
     
      }
  //end of genParticle loop
      if (GEN_MUON_LIST.size()>1){
        sort (GEN_MUON_LIST.begin(), GEN_MUON_LIST.end(), PtComparator<GenParticle>());
        const GenParticle& mua1 = GEN_MUON_LIST[0];
        const GenParticle& mua2 = GEN_MUON_LIST[1];
        TLorentzVector x1, x2;
        x1.SetPtEtaPhiE (mua1.pt, mua1.eta, mua1.phi, mua1.energy);
        x2.SetPtEtaPhiE (mua2.pt, mua2.eta, mua2.phi, mua2.energy);
        TLorentzVector sum = x1+x2 ;
        double MASS = sum.M();
        AnaUtil::fillHist1D ("genInvMass", MASS, 1.0);
      }
      
      //AnaUtil::fillHist1D ("muon_count", genMuonList.size(), 1.0);  
      if (genMuonList.size() > 1) 
        sort(genMuonList.begin(), genMuonList.end(), PtComparator<GenParticle>());
      //      if (genMuonList.size() > 0) {
      //   double pt = genMuonList[0].pt;
      //  AnaUtil::fillHist1D ("genMuon_pt1", pt, 1.0);
      // }
      if (genTauList.size() > 1)
        sort(genTauList.begin(), genTauList.end(), PtComparator<GenParticle>());
      if (genMuoaList.size() > 1)
        sort(genMuoaList.begin(), genMuoaList.end(), PtComparator<GenParticle>());
      if (genMuobList.size() > 1)
        sort(genMuobList.begin(), genMuobList.end(), PtComparator<GenParticle>());
      if (genHList.size() > 1)
        sort(genHList.begin(), genHList.end(), PtComparator<GenParticle>());
      if (genWList.size() > 1)
        sort(genWList.begin(), genWList.end(), PtComparator<GenParticle>());
 
     if (genMuonList.size()>1 && genTauList.size()>0){         
         const GenParticle& muona=genMuonList[0];
         const GenParticle& muonb=genMuonList[1];
         const GenParticle& taua =genTauList[0];
         if (muona.pt<15) return;
         if (muonb.pt<10) return;
         TLorentzVector v1,v2,t1;
         v1.SetPtEtaPhiE (muona.pt, muona.eta, muona.phi, muona.energy);
         v2.SetPtEtaPhiE (muonb.pt, muonb.eta, muonb.phi, muonb.energy);
         t1.SetPtEtaPhiE (taua.pt , taua.eta , taua.phi , taua.energy );
         TLorentzVector z= v1+v2;
         double mass= z.M();
         if (mass < 20) return;
         double mva_gen= z.Pt()/(v1.Pt()+v2.Pt());
  
         AnaUtil::fillHist1D("mva_gen", mva_gen, 1.0);
         AnaUtil::fillHist2D("muon_correlation", v1.Pt(), v2.Pt(), 1.0);
	 AnaUtil::fillHist1D("oddInvMass",mass, 1.0);
      }
      if (genMuoaList.size()>0 && genMuobList.size()>0 && genTauList.size()>0) {
      
       	TLorentzVector M1,M2,W1,H1,t1;
        M1.SetPtEtaPhiE (genMuoaList[0].pt, genMuoaList[0].eta, genMuoaList[0].phi, genMuoaList[0].energy);
        M2.SetPtEtaPhiE (genMuobList[0].pt, genMuobList[0].eta, genMuobList[0].phi, genMuobList[0].energy);
        W1.SetPtEtaPhiE (genWList[0].pt, genWList[0].eta, genWList[0].phi, genWList[0].energy);
        H1.SetPtEtaPhiE (genHList[0].pt, genHList[0].eta, genHList[0].phi, genHList[0].energy);
        t1.SetPtEtaPhiE (genTauList[0].pt, genTauList[0].eta, genTauList[0].phi, genTauList[0].energy);
        AnaUtil::fillHist1D("WH_DR", AnaUtil::deltaR(W1, H1), 1.0);
        AnaUtil::fillHist1D("Ht1_DR",AnaUtil::deltaR(t1, H1), 1.0);
        AnaUtil::fillHist1D("HM2_DR",AnaUtil::deltaR(M2, H1), 1.0);
        AnaUtil::fillHist1D("WM1_DR",AnaUtil::deltaR(M1, W1), 1.0);
        AnaUtil::fillHist1D("t1M2_DR",AnaUtil::deltaR(t1,M2), 1.0);
        AnaUtil::fillProfile("WH_Wpt",W1.Pt(),AnaUtil::deltaR(W1, H1), 1.0);
        AnaUtil::fillProfile("t1M2_Hpt",H1.Pt(),AnaUtil::deltaR(t1,M2), 1.0);
        AnaUtil::fillProfile("WM1_Wpt", W1.Pt(),AnaUtil::deltaR(W1,M1), 1.0);
        AnaUtil::fillHist1D("DPhi_WH",AnaUtil::deltaPhi(W1,H1), 1.0);
        AnaUtil::fillHist1D("DEta_WH",W1.Eta()-H1.Eta(), 1.0);
        AnaUtil::fillHist1D("DPhi_t1M2",AnaUtil::deltaPhi(t1,M2), 1.0);
        AnaUtil::fillHist1D("DEta_t1M2",t1.Eta()-M2.Eta(), 1.0);
      }

      
    }
  //end of generator level   
    /* if (genTau_H_List.size()>0 && genMuon_Z_List.size()>1){
      signal_ZH = signal_ZH + 1;
    }

    if (genMuonList.size()>1 && genTauList.size()>0){
      signal= signal+1;
    }*/
}
void MuTauTau::computeDeltaR(const vector<Muon>& muoList, const vector<Tau>& tauList ){
 if(muoList.size() > 0 && tauList.size() > 1){
  const Muon& muoa = muoList[0];

  if (muoa.pt < 20) return;

  TLorentzVector m1;
  m1.SetPtEtaPhiE(muoa.pt, muoa.eta, muoa.phi, muoa.energy); 

 }
}
/*void AnaBase::GenDetLevelMatch (const vector<Muon>& muoList,
                                const vector<GenParticle>& genMuonList){
for (int i = 0; i < 2; i++){
    double minDR = 990;
    int index = -1;
    for(int j = 0; j < genMuonList.size();j++){
       double deltaDR = sqrt(pow((genMuonList[j].eta - muolist[i].eta),2) + pow((genMuonList[j].phi - muolist[i].phi),2));
       if(deltaDR<minDr){
          minDR = deltaDR;
          index = j; //gen label index
          }
       }
    if(index < 0) continue;
    
    if(i = 0){
       AnaUtil::fillHist1D("deltaPt1",muolist[0].pt - genMuonList[index].pt);
       AnaUtil::fillHist1D("deltaPt1",muolist[0].eta - genMuonList[index].eta);
       AnaUtil::fillHist1D("deltaPt1",muolist[0].phi - genMuonList[index].phi);
    }
    else{
       AnaUtil::fillHist1D("deltaDR1",minDR);
       AnaUtil::fillHist1D("deltaPt1",muolist[1].pt - genMuonList[index].pt);
       AnaUtil::fillHist1D("deltaPt1",muolist[1].eta - genMuonList[index].eta);
       AnaUtil::fillHist1D("deltaPt1",muolist[1].phi - genMuonList[index].phi);
       }  
    }
}
*/
void MuTauTau::selectEvent()
{
  int nsmuon = muoList.size();
  int nstau = tauList.size();

  // require at least 1 muon
  if (nsmuon < 1) return;
  nEvtSel[3] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 3, _puevWt);

  // muon Pt
  const Muon& muoa = muoList.at(0);
  if (muoa.pt <= AnaUtil::cutValue(_evselCutMap, "ptMuon")) return;
  nEvtSel[4] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 4, _puevWt);

  AnaUtil::fillHist1D("muonPt1_stage1",  muoa.pt, _puevWt);
  AnaUtil::fillHist1D("muonEta1_stage1", muoa.eta, _puevWt);

  // 1 hadronic tau
  if (nstau < 1) return;
  nEvtSel[5] += _puevWt;;
  AnaUtil::fillHist1D("ecounter", 5, _puevWt);

  const Tau& taua = tauList.at(0);

  // Tau Pt
  if (taua.pt <= AnaUtil::cutValue(_evselCutMap, "ptTau1")) return;
  nEvtSel[6] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 6, _puevWt);

  AnaUtil::fillHist1D("muonPt1_stage2", muoa.pt, _puevWt);
  AnaUtil::fillHist1D("muonEta1_stage2", muoa.eta, _puevWt);
  AnaUtil::fillHist1D("tauPt1_stage1", taua.pt, _puevWt);
  AnaUtil::fillHist1D("tauEta1_stage1", taua.eta, _puevWt);

  // 2 hadronic tau
  if (nstau < 2) return;
  nEvtSel[7] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 7, _puevWt);

  const Tau& taub = tauList.at(1);
  if (taub.pt <= AnaUtil::cutValue(_evselCutMap, "ptTau2")) return;
  nEvtSel[8] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 8, _puevWt);

  AnaUtil::fillHist1D("muonPt1_stage3", muoa.pt, _puevWt);
  AnaUtil::fillHist1D("muonEta1_stage3", muoa.eta, _puevWt);
  AnaUtil::fillHist1D("tauPt1_stage2", taua.pt, _puevWt);
  AnaUtil::fillHist1D("tauEta1_stage2", taua.eta, _puevWt);
  AnaUtil::fillHist1D("tauPt2_stage1", taub.pt, _puevWt);
  AnaUtil::fillHist1D("tauEta2_stage1", taub.eta, _puevWt);

  // no b-tagged jet
  if (bjetList.size() > 0) return;
  nEvtSel[9] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 9, _puevWt);

  // Opposite charge for tau_1 & tau_2
  if ( (taua.charge + taub.charge) != 0) return;
  nEvtSel[10] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 10, _puevWt);
  
  // no good electron
  if (eleList.size() > 0) return;
  nEvtSel[11] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 11, _puevWt);

  TLorentzVector M1, T1, T2 ;
  M1.SetPtEtaPhiE(muoa.pt, muoa.eta, muoa.phi, muoa.energy);
  T1.SetPtEtaPhiE(taua.pt, taua.eta, taua.phi, taua.energy);
  T2.SetPtEtaPhiE(taub.pt, taub.eta, taub.phi, taub.energy);  
  double drM1T1 = AnaUtil::deltaR(M1, T1);
  double drM1T2 = AnaUtil::deltaR(M1, T2);
  double drT1T2 = AnaUtil::deltaR(T1, T2);
  AnaUtil::fillHist1D("deltaR_mutau1", drM1T1, _puevWt);
  AnaUtil::fillHist1D("deltaR_mutau2", drM1T2, _puevWt);
  AnaUtil::fillHist1D("deltaR_tau1tau2", drT1T2, _puevWt);
  
  if (drM1T1 <= AnaUtil::cutValue(_evselCutMap, "drMuTau1")) return;
  nEvtSel[12] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 12, _puevWt);

  if (drM1T2 <= AnaUtil::cutValue(_evselCutMap, "drMuTau2")) return;
  nEvtSel[13] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 13, _puevWt);
 
  //temp entry, to be removed
  //if (muoList.size() > 1 && (muo.charge + muoList[1].charge == 0)) return;

  if (drT1T2 <= AnaUtil::cutValue(_evselCutMap, "drTau1Tau2")) return;
  nEvtSel[14] += _puevWt;
  AnaUtil::fillHist1D("ecounter", 14, _puevWt);

  AnaUtil::fillHist1D("nMuon", nsmuon, _puevWt);
  AnaUtil::fillHist1D("muonPt1_allcut", muoa.pt, _puevWt);
  AnaUtil::fillHist1D("muonEta1_allcut", muoa.eta, _puevWt);
  if (nsmuon > 1) {
    const Muon& muob = muoList.at(1);    
    AnaUtil::fillHist1D("muonPt2", muob.pt, _puevWt);
    TLorentzVector m1, m2;
    m1.SetPtEtaPhiE(muoa.pt, muoa.eta, muoa.phi, muoa.energy);
    m2.SetPtEtaPhiE(muob.pt, muob.eta, muob.phi, muob.energy);
    TLorentzVector z = m1 + m2;
    double mass = z.M();
    AnaUtil::fillHist1D("mumuInv", mass, _puevWt);
  }
  if (nstau > 2) AnaUtil::fillHist1D("tauPt3", tauList.at(2).pt, _puevWt);
  AnaUtil::fillHist1D("nTau", nstau, _puevWt);
  int n5tau = 0;
  for (unsigned int i = 0; i < tauList.size(); ++i) {
    if (tauList.at(i).pt >= 5) n5tau = n5tau + 1;
  }
  AnaUtil::fillHist1D("nTau_gt5", n5tau, _puevWt);
  AnaUtil::fillHist1D("tauPt1_allcut", taua.pt, _puevWt);
  AnaUtil::fillHist1D("tauEta1_allcut", taua.eta, _puevWt);
  AnaUtil::fillHist1D("tauPt2_allcut", taub.pt, _puevWt);
  AnaUtil::fillHist1D("tauEta2_allcut", taub.eta, _puevWt);

  TLorentzVector t = T1 + T2;
  double invmass = t.M();
  AnaUtil::fillHist1D("tautauInvmass", invmass, _puevWt); 

  const MET* mt = dynamic_cast<MET*>(metA->At(0));
  assert(mt);
  AnaUtil::fillHist1D("met", mt->met, _puevWt);
  double dphi = AnaUtil::deltaPhi(muoa.phi, mt->metphi);
  AnaUtil::fillHist1D("diffPhi_Mu_Met", dphi, _puevWt);
  AnaUtil::fillHist1D("diffPhi_Mu_Tau1", drM1T1, _puevWt);
  AnaUtil::fillHist1D("DiTauPt", t.Pt(), _puevWt);
  AnaUtil::fillHist1D("diffPhi_Mu_DiTau", AnaUtil::deltaPhi(M1, t), _puevWt);
  double var = t.Pt()/(T1.Pt() + T2.Pt());
  AnaUtil::fillHist1D("var1", var, _puevWt); 
  double mass1 = sqrt(pow(M1.Pt() + mt->met,2) - pow(M1.Pt()*cos(M1.Phi()) + mt->met * cos(mt->metphi),2) - 
		      pow(M1.Pt()*sin(M1.Phi()) + mt->met* sin(mt->metphi),2));   
  double mass2 = sqrt(pow(T1.Pt() + mt->met,2) - pow(T1.Pt()*cos(T1.Phi()) + mt->met * cos(mt->metphi),2) - 
		      pow(T1.Pt()*sin(T1.Phi()) + mt->met* sin(mt->metphi),2));   
  double mass3 = sqrt(pow(T2.Pt() + mt->met,2) - pow(T2.Pt()*cos(T2.Phi()) + mt->met * cos(mt->metphi),2) - 
		      pow(T2.Pt()*sin(T2.Phi()) + mt->met* sin(mt->metphi),2));   
  AnaUtil::fillHist1D("MT_Mu_Met", mass1, _puevWt);
  AnaUtil::fillHist1D("MT_Tau1_Met", mass2, _puevWt);
  AnaUtil::fillHist1D("MT_Tau2_Met", mass3, _puevWt);

  // Dump seleted event information
  // Event and object information
  vector<map<string, double> > maps;
  maps.push_back(_vtxCutMap);
  maps.push_back(_tauCutMap);
  maps.push_back(_muonCutMap);
  maps.push_back(_electronCutMap);
  maps.push_back(_bjetCutMap);
  dumpEvent(maps, _evLog);

  // event selection information
  _evLog << setprecision(3);
  _evLog << "nVertex = " << vtxList.size() << endl
         << "nMuon = " << nsmuon << endl
         << "pT Muon 1 = " << setw(8) << muoa.pt << " GeV" << endl;
  if (nsmuon > 1) _evLog << "pT Muon 2 = " << setw(8) << muoList.at(1).pt << " GeV" << endl;
  _evLog << "nTau = " << nstau << endl
         << "pT Tau 1 = " << setw(8) << taua.pt << " GeV" << endl
         << "pT Tau 2 = " << setw(8) << taub.pt << " GeV" << endl
         << "DR(m,t1) = " << setw(7) << drM1T1 << endl
         << "DR(m,t2) = " << setw(7) << drM1T2 << endl
         << "DR(t1,t2) = " << setw(7) << drT1T2 << endl
         << "mass(t1,t2) = " << setw(8) << invmass << " GeV" << endl
         << "MET = " << setw(8) << mt->met << " GeV"
         << endl;
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void MuTauTau::endJob() 
{  
  _fLog << setprecision(1);
  _fLog << "============================================================" << endl
        << "Statistics for W->munu, H->tau+tau-; tau->hadron,tau->hadron" << endl
        << "============================================================" << endl  
        << "       Total events processed:" << setw(10) << nEvtSel[0] << endl 
        << "     after Trigger selectionD:" << setw(10) << nEvtSel[1] << endl        
        << "       >= 1 Good event vertex:" << setw(10) << nEvtSel[2] << endl        
        << "          >= 1 selected muons:" << setw(10) << nEvtSel[3] << endl
        << "           Muon Pt > " << AnaUtil::cutValue(_evselCutMap, "ptMuon") << " GeV:"
                                            << setw(10) << nEvtSel[4] << endl
        << "           >= 1 selected taus:" << setw(10) << nEvtSel[5] << endl
        << "    Highest Tau Pt > " << AnaUtil::cutValue(_evselCutMap, "ptTau1") << " GeV:"
                                            << setw(10) << nEvtSel[6] << endl
	<< "           >= 2 selected taus:" << setw(10) << nEvtSel[7] << endl
        << "        2nd Tau Pt > " << AnaUtil::cutValue(_evselCutMap, "ptTau2") << " GeV:"
                                            << setw(10) << nEvtSel[8] << endl
        << "              no tagged b-jet:" << setw(10) << nEvtSel[9] << endl
        << "          tau+tau charge == 0:" << setw(10) << nEvtSel[10] << endl
        << "             no good electron:" << setw(10) << nEvtSel[11] << endl
        << "      Muon Tau1 overlap > " << AnaUtil::cutValue(_evselCutMap, "drMuTau1") << ":"
                                            << setw(10) << nEvtSel[12] << endl
        << "      Muon Tau2 overlap > " << AnaUtil::cutValue(_evselCutMap, "drMuTau2") << ":"
                                            << setw(10) << nEvtSel[13] << endl
        << "      Tau1 Tau2 overlap > " << AnaUtil::cutValue(_evselCutMap, "drTau1Tau2") << ":"
                                            << setw(10) << nEvtSel[14] << endl;
  _fLog << resetiosflags(ios::fixed);

  closeFiles();

  _histf->cd();
  _histf->Write();
  _histf->Close();
  delete _histf;
}
// -------------------------------------------------------------------------------
// Poor man's way of a datacard. Each line between the 'START' and 'END' tags
// is read in turn, split into words, where the first element is the 'key' and
// the rest the value(s). If more than one values are present they are 
// stored in a vector. No safety mechanism is in place. Any line with an unknown 
// key is skipped. Comments lines should start with either '#' or '//', preferably
// in the first column. Empty lines are skipped. The file containing the datacards 
// is passed as the only argument of the program, there is no default datacard
// -------------------------------------------------------------------------------
bool MuTauTau::readJob(const string& jobFile, int& nFiles)
{
  AnaBase::readJob(jobFile, nFiles);

  static const int BUF_SIZE = 256;

  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }

  // note that you must use a pointer (reference!) to the cut map
  // in order to avoid scope related issues
  map<string, map<string, double>* > hmap;
  hmap.insert(pair<string, map<string, double>* >("vtxCutList", &_vtxCutMap));
  hmap.insert(pair<string, map<string, double>* >("electronCutList", &_electronCutMap));
  hmap.insert(pair<string, map<string, double>* >("muonCutList", &_muonCutMap));
  hmap.insert(pair<string, map<string, double>* >("tauCutList", &_tauCutMap));
  hmap.insert(pair<string, map<string, double>* >("bjetCutList", &_bjetCutMap));
  hmap.insert(pair<string, map<string, double>* >("evselCutList", &_evselCutMap));

  char buf[BUF_SIZE];
  vector<string> tokens;
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   

    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words
    AnaUtil::tokenize(line, tokens);
    int size = tokens.size();
    string key = tokens[0];

    map<string, map<string, double>* >::const_iterator pos = hmap.find(key);
    if (pos != hmap.end()) {
      map<string, double>* m = pos->second;        
      m->clear();
      for (int i = 1; i < size; ++i) {
        vector<string> cutstr;
        // Split the line into words
	AnaUtil::tokenize(tokens[i], cutstr, "=");
        m->insert( pair<string,double>(cutstr[0], atof(cutstr[1].c_str())));
      }
    }    
    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();

  return true;
}
void MuTauTau::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("vtxCutList", _vtxCutMap));
  hmap.insert(pair<string, map<string, double> >("electronCutList", _electronCutMap));
  hmap.insert(pair<string, map<string, double> >("muonCutList", _muonCutMap));
  hmap.insert(pair<string, map<string, double> >("tauCutList", _tauCutMap));
  hmap.insert(pair<string, map<string, double> >("bjetCutList", _bjetCutMap));
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));

  for (map<string, map<string, double> >::const_iterator it  = hmap.begin(); 
                                                         it != hmap.end(); ++it)  
  {
    os << ">>> " << it->first << endl; 
    map<string, double> m = it->second;
    os << setprecision(2);
    for (map<string,double>::const_iterator jt  = m.begin(); 
                                            jt != m.end(); ++jt)  
      os << setw(16) << jt->first << ": " 
         << setw(7) << jt->second << endl;
    os << endl; 
  }
}
