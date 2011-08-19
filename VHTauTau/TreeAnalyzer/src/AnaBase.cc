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

#include "AnaBase.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

class VertexComparator {
public:
  bool operator()(const Vertex &a, const Vertex &b) const {
    return a.sumPt > b.sumPt;
  }
};

template <class T>
class PtComparator {
public:
  bool operator()(const T &a, const T &b) const {
    return a.pt > b.pt;
  }
};
// -----------
// Constructor
// -----------
AnaBase::AnaBase(const string& filename)
  : _chain(new TChain("treeCreator/vhtree")),
    _histf(0),
          eventA(new TClonesArray("Event")),
        triggerA(new TClonesArray("Trigger")),
         vertexA(new TClonesArray("Vertex")),    
            tauA(new TClonesArray("Tau")),    
       electronA(new TClonesArray("Electron")),
           muonA(new TClonesArray("Muon")),    
            jetA(new TClonesArray("Jet")),    
            metA(new TClonesArray("MET")),    
    genParticleA(new TClonesArray("GenParticle")),  
         genJetA(new TClonesArray("GenJet")),   
         genMetA(new TClonesArray("GenMET")),
       //genEventA(new TClonesArray("GenEvent")),   
      _logOption(0),
         _maxEvt(0)
{
  cout << setiosflags(ios::fixed); 
  cout << "=== Start of Analysis === " << endl;
}
// ----------
// Destructor
// ----------
AnaBase::~AnaBase() 
{
  clearEvent();

  delete eventA;
  delete triggerA;
  delete vertexA;
//  delete genEventA;
  delete tauA;
  delete electronA;
  delete muonA;
  delete jetA;
  delete metA;
  delete genParticleA;
  delete genJetA;
  delete genMetA;
}
// ------------------------
// Clear the clones arrays
// ------------------------
void AnaBase::clearEvent() 
{
  if (eventA) eventA->Clear();
//  if (genEventA) genEventA->Clear();
  if (triggerA) triggerA->Clear();
  if (vertexA) vertexA->Clear();
  if (tauA) tauA->Clear();
  if (electronA) electronA->Clear();
  if (muonA) muonA->Clear();
  if (jetA) jetA->Clear();
  if (metA) metA->Clear();
  if (genParticleA) genParticleA->Clear();
  if (genJetA) genJetA->Clear();
  if (genMetA) genMetA->Clear();


  n_vertex      = 0;  
  n_tau         = 0;
  n_electron    = 0;
  n_muon        = 0;
  n_jet         = 0;
  n_met         = 0;
  n_genparticle = 0;
  n_genjet      = 0;
  n_genmet      = 0;
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool AnaBase::beginJob() 
{
  for (unsigned int i = 0; i < NEL(nEvt); i++) {
    nEvt[i] = 0;
  }
  // Open the output ROOT file
  _histf = TFile::Open(_histFile.c_str(), "RECREATE");
  bookHistograms();

  setAddresses();
  //enableBranches();
  nEvents = static_cast<int>(_chain->GetEntries()); 
  if (nEvents <= 0) {
    cerr << "******* nEvents = " << nEvents << ", returning!" << endl;
    return false;
  }
  if (_maxEvt > 0) nEvents = _maxEvt;
  cout << " ===== # of events to analyse, nEvents = " << nEvents << endl;

  openFiles();

  return true;
}
// ---------------
// Book histograms
// ---------------
void AnaBase::bookHistograms() 
{
  //static const double pi = TMath::Pi();
  new TH1F("muonPt", "Muon Pt Distribution", 100, -0.5, 199.5);
  new TH1F("elecPt", "Electron Pt Distribution", 100, -0.5, 199.5);
}
// -------------------
// The main event loop
// -------------------
void AnaBase::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;

  int nPrint = std::max(10000, nEvents/1000);

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  for (int ev = 0; ev < nEvents; ev++) {
    bool select = false;

    clearEvent();
    int lflag = _chain->LoadTree(ev); 
    int nentries = getEntry(lflag);    // returns bytes read
    ++nEvt[0];
    
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
           << " Events proc. " << setw(8) << ev << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
       cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
            << " ==> " << _chain->GetCurrentFile()->GetName() 
            << " <<< Run# " << run 
            << " Event# " << setw(8) << event << " >>> " 
            << " Events proc. " << setw(8) << ev << endl;

    if (_logOption >> 0 & 0x1) 
    _fLog << "n_tau: "<< n_tau
          << ", n_muon: "<< n_muon
          << ", n_jet: " << n_jet
          << ", n_vertex: " << n_vertex
          << ", n_met: " << n_met
          << ", n_electron: " << n_electron << endl;

    vector<Vertex> vtxList;
    vector<Muon> muoList;
    vector<Electron> eleList;
    vector<Tau> tauList;
    vector<Jet> bjetList;

    // Let's look at the event vertex  
    _fLog << setprecision(3);
    if (n_vertex && (_logOption >> 1 & 0x1)) {
      _fLog << "=>> Vertices: " << n_vertex << endl
            << "indx     ndf     dxy       z   sumPt    chi2   ntrks  ntrkw05     sbit" << endl; 
    }
    for (int indx = 0; indx < n_vertex; ++indx) {
      const Vertex* vtx = dynamic_cast<Vertex*>(vertexA->At(indx));
      if (!vtx) continue;

      double dxy = std::sqrt(pow(vtx->x, 2) + pow(vtx->y, 2));
      int sbit = 0;
      if (vtx->ndf <= _vtxCutMap["ndf"])       sbit |= (1 << 0);
      if (dxy >= _vtxCutMap["dxy"])            sbit |= (1 << 1);
      if (std::abs(vtx->z) >= _vtxCutMap["z"]) sbit |= (1 << 2);

      if (_logOption >> 1 & 0x1) {
        _fLog << setw(4) << indx
              << setw(8) << vtx->ndf
              << setw(8) << dxy
              << setw(8) << vtx->z 
              << setw(8) << vtx->sumPt
              << setw(8) << vtx->chi2
              << setw(8) << vtx->ntracks
              << setw(9) << vtx->ntracksw05;
        AnaUtil::bit_print(sbit, 8, _fLog);
      }
      if (sbit) continue;
      vtxList.push_back(*vtx);
    }
    if (vtxList.size() > 1) 
      std::sort(vtxList.begin(), vtxList.end(), VertexComparator());
#if 0
    if (n_genparticle) {
      for (int indx = 0; indx < n_genparticle; ++indx) {
        const GenParticle* gp = dynamic_cast<GenParticle*>(genParticleA->At(indx));
        if (!gp) continue;
        _fLog << "pdgId: " << gp->pdgId << endl;

	vector<int> m = gp->motherIndices;
	_fLog << "# of mother: "  << m.size() << endl;
        for (size_t i = 0; i < m.size(); ++i) {
          int mi = m[i];
          _fLog << "index: " << mi << endl;
          if (mi >= n_genparticle) continue;
          const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(mi));
          _fLog << "\tpdgId: " << mgp->pdgId << endl; 
        }

	vector<int> d = gp->daughtIndices;
	_fLog << "# of daughter: " << d.size() << endl;
        for (size_t i = 0; i < d.size(); ++i) {
	  int di = d[i];
          _fLog << "index: " << di << endl;
          if (di >= n_genparticle) continue;
          const GenParticle* dgp = dynamic_cast<GenParticle*>(genParticleA->At(di));
          _fLog << "\tpdgId: " << dgp->pdgId << endl; 
        }
      }
    }
#endif
    // Let's look at the Tau collection
    if (n_tau && (_logOption >> 2 & 0x1)) {
      _fLog << "=>> Taus: " << n_tau << endl
            << "indx     Eta     Phi      Pt  Energy"
  	    << "   lchPt   lnpPt    lpPt DMF  LI aMT aEL     sbit" 
	    << endl; 
    }
    for (int indx = 0; indx < n_tau; ++indx) {
      const Tau* tau = dynamic_cast<Tau*>(tauA->At(indx));
      if (!tau) continue;

      int sbit = 0;
      if (abs(tau->eta) >= _tauCutMap["eta"]) sbit |= (1 << 0); 
      if (tau->decayModeFinding     <= 0.5)   sbit |= (1 << 1); 
      if (tau->looseIsolation       <= 0.5)   sbit |= (1 << 2); 
      if (tau->againstMuonTight     <= 0.5)   sbit |= (1 << 3); 
      if (tau->againstElectronLoose <= 0.5)   sbit |= (1 << 4); 

      if (_logOption >> 2 & 0x1) {
        _fLog << setprecision(2)
              << setw(4) << indx 
              << setw(8) << tau->eta
              << setw(8) << tau->phi
              << setw(8) << tau->pt
              << setw(8) << tau->energy
              << setw(8) << tau->leadChargedParticlePt
              << setw(8) << tau->leadNeutralParticlePt
              << setw(8) << tau->leadParticlePt
              << setprecision(1)
              << setw(4) << tau->decayModeFinding
              << setw(4) << tau->looseIsolation
              << setw(4) << tau->againstMuonTight
              << setw(4) << tau->againstElectronLoose;
        AnaUtil::bit_print(sbit, 8, _fLog);
      }
      if (sbit) continue;
      tauList.push_back(*tau);
    }
    if (tauList.size() > 1) 
      std::sort(tauList.begin(), tauList.end(), PtComparator<Tau>());
   
    // Let's look at the Muon collection
    if (n_muon && (_logOption >> 3 & 0x1)) {
      _fLog << "=>> Muons: " << n_muon << endl
            << "indx     Eta     Phi      Pt       P"
  	    << "      D0   D0Err      Dz   DzErr" 
            << "  relIso pixHits trkHits   gChi2      dB     sbit"
	    << endl; 
    }
    for (int indx = 0; indx < n_muon; ++indx) {
      const Muon* muon = dynamic_cast<Muon*>(muonA->At(indx));
      if (!muon) continue;

      int sbit = 0;
      if (!muon->isTrackerMuon)                                     sbit |= (1 << 0); 
      if (std::abs(muon->eta) >= _muonCutMap["eta"])                sbit |= (1 << 1);
      if (muon->relIso >= _muonCutMap["relIso"])                    sbit |= (1 << 2);
      if ((muon->pixHits + muon->trkHits) <= _muonCutMap["ptHits"]) sbit |= (1 << 3);
      if (muon->globalChi2 >= _muonCutMap["globalChi2"])            sbit |= (1 << 4);
      if (std::abs(muon->trkD0) >= _muonCutMap["trkD0"])            sbit |= (1 << 5);
      if (std::abs(muon->dB) >= _muonCutMap["dB"])                  sbit |= (1 << 6);

      if (_logOption >> 3 & 0x1) {
        _fLog << setw(4) << indx 
              << setw(8) << muon->eta
              << setw(8) << muon->phi
              << setw(8) << muon->pt
              << setw(8) << muon->p
              << setw(8) << muon->trkD0
              << setw(8) << muon->trkD0Error
              << setw(8) << muon->trkDz
              << setw(8) << muon->trkDzError
              << setw(8) << muon->relIso
              << setw(8) << muon->pixHits
              << setw(8) << muon->trkHits
              << setw(8) << muon->globalChi2
              << setw(8) << muon->dB;
        AnaUtil::bit_print(sbit, 8, _fLog);
      }

      // Now apply cuts
      if (sbit) continue;
      muoList.push_back(*muon);
      fillHist1D("muonPt", muon->pt, 1.0);
    }
    if (muoList.size() > 1) 
      std::sort(muoList.begin(), muoList.end(), PtComparator<Muon>());

    // Let's look at the Electron collection now
    if (n_electron && (_logOption >> 4 & 0x1)) {
      _fLog << "=>> Electrons: " << n_electron << endl
            << "indx     Eta     Phi      Pt  Energy"
  	    << "   scEta   dB      eleId     sbit" 
	    << endl; 
    } 
    for (int indx = 0; indx < n_electron; ++indx) {
      const Electron* elec = dynamic_cast<Electron*>(electronA->At(indx));
      if (!elec) continue;

      int sbit = 0;
      if (elec->pt <= _electronCutMap["pt"])                    sbit |= (1 << 0);
      if (std::abs(elec->eta) >= _electronCutMap["eta"])        sbit |= (1 << 1);
      if (elec->simpleEleId95cIso <= _electronCutMap["eleId"])  sbit |= (1 << 2);
      if (!elec->hasGsfTrack)                                   sbit |= (1 << 3);
      if (std::abs(elec->dB) >= _electronCutMap["dB"])          sbit |= (1 << 4);
      if ( std::abs(elec->scEta) >= _electronCutMap["scEtaLow"] 
        && std::abs(elec->scEta) <= _electronCutMap["scEtaUp"]) sbit |= (1 << 5);

      if (_logOption >> 4 & 0x1) {
        _fLog << setw(4) << indx 
              << setw(8) << elec->eta
              << setw(8) << elec->phi
              << setw(8) << elec->pt
              << setw(8) << elec->energy
              << setw(8) << elec->scEta
              << setw(8) << elec->dB
              << setw(8) << elec->simpleEleId95cIso;
        AnaUtil::bit_print(sbit, 8, _fLog);
      }
      // Now apply cuts
      if (sbit) continue;
      eleList.push_back(*elec);
      fillHist1D("elecPt", elec->pt, 1.0);
    }
    if (eleList.size() > 1) 
      std::sort(eleList.begin(), eleList.end(), PtComparator<Electron>());

    // Let's look at the Jet collection
    if (n_jet && (_logOption >> 5 & 0x1)) {
      _fLog << "=>> Jets: " << n_jet << endl
            << "indx     Eta     Phi      Pt  Energy"
            << "    TCHE    TCHP     sbit"
	    << endl; 
    } 
    for (int indx = 0; indx < n_jet; ++indx) {
      const Jet* jt = dynamic_cast<Jet*>(jetA->At(indx));
      if (!jt) continue;

      int sbit = 0;
      if (jt->eta >= _bjetCutMap["eta"])                             sbit |= (1 << 0);
      if (jt->pt <= _bjetCutMap["pt"])                               sbit |= (1 << 1);
      if (jt->trackCountingHighEffBTag <= _bjetCutMap["trackCount"]) sbit |= (1 << 2);

      if (_logOption >> 5 & 0x1) {
        _fLog << setw(4) << indx 
              << setw(8) << jt->eta
              << setw(8) << jt->phi
              << setw(8) << jt->pt
              << setw(8) << jt->energy
              << setprecision(1)
              << setw(8) << jt->trackCountingHighEffBTag  
              << setw(8) << jt->trackCountingHighPurBTag
              << setw(8) << jt->jetProbabilityBTag
              << setw(8) << jt->jetBProbabilityBTag;
        AnaUtil::bit_print(sbit, 8, _fLog);
      }
      // Now apply cuts
      if (sbit) continue;
      bjetList.push_back(*jt);
    }
    if (bjetList.size() > 1) 
      std::sort(bjetList.begin(), bjetList.end(), PtComparator<Jet>());

    if (_logOption)
    _fLog << "n_vertex_good: "     << vtxList.size()
          << ", n_muon_selected: " << muoList.size()
          << ", n_electron: "      << eleList.size()
          << ", n_tau_selected: "  << tauList.size()
          << ", n_bjet: "          << bjetList.size()
          << endl;

    // presence of > 1 good vertex
    if (vtxList.size() < 1) continue;
    ++nEvt[1];

    // 2 muons
    if (muoList.size() < 2) continue;
    ++nEvt[2];

    // muon Pt
    const Muon& muoa = muoList[0];
    const Muon& muob = muoList[1];
    if (muoa.pt <= 15) continue;
    ++nEvt[3];

    if (muob.pt <= 10) continue;
    ++nEvt[4];

    // 1 hadronic tau
    if (tauList.size() < 1) continue;
    ++nEvt[5];

    const Tau& tau = tauList[0];

    // Tau Pt
    if (tau.pt <= 15) continue;
    ++nEvt[6];

    // no b-tagged jet
    if (bjetList.size()) continue;
    ++nEvt[7];

    // Opposite charge for mu_2 & tau
    if ( (muob.charge + tau.charge) != 0) continue;
    ++nEvt[8];

    // Same charge for muon_1 & muon_2 (to remove DY)
    if ( (muoa.charge + muob.charge) == 0) continue;
    ++nEvt[9];
  }  
  // Analysis is over
  endJob();
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void AnaBase::endJob() 
{
  
  const string tags[] = 
  {
    "Total events processed",
    ">= 1 Good event vertex",
    ">= 2 selected muons",
    "First Muon Pt> 15 GeV cut",
    "Second Muon Pt>10 GeV cut",
    ">= 1 selected tau",
    "First Tau Pt>15 GeV cut",
    "no tagged b-jet",
    "mu+tau charge == 0",
    "mu+mu charge != 0"
  };
  for (unsigned int i = 0; i < NEL(nEvt); i++) {
    _fLog << setw(64) << tags[i] << setw(8) << nEvt[i] << endl;
  }
  _fLog << resetiosflags(ios::fixed);

  closeFiles();

  _histf->cd();
  _histf->Write();
  _histf->Close();
  delete _histf;
}

// ----------------------------------------------------------
// Perform event selection, For selection of Z -> e+e- events
// we need,
//   - > 0 Tight electron
//   - event within e+e- invariant mass window
// ----------------------------------------------------------
bool AnaBase::selectEvent() 
{
  return true;
}
// ------------------------------------
// Get Run number for the present event
// ------------------------------------
int AnaBase::getRunNumber() const 
{
  const Event* event = dynamic_cast<Event*>(eventA->At(0));
  assert(event);
  return event->run;
}    
// ---------------------------------
// Add input Root files to the chain
// ---------------------------------
void AnaBase::setInputFile(const string& fname) 
{
  size_t found = fname.find("root:");
  if (found == string::npos && gSystem->AccessPathName(fname.c_str())) {
     cerr << "=>> Warning: File <<" << fname << ">> was not found!!" << endl;
     return;  
  }
  _chain->AddFile(fname.c_str(), -1);
}
// ---------------------------------------
// Get total number of events in the chain
// --------------------------------------
int AnaBase::getEntries() const 
{
  return static_cast<int>(_chain->GetEntries());
}
// ------------------------------------------------------------------------
// Convenience routine for filling 1D histograms. We rely on root to keep 
// track of all the histograms that are booked all over so that we do not 
// have to use any global variables to save the histogram pointers. Instead, 
// we use the name of the histograms and gROOT to retrieve them from the 
// Root object pool whenever necessary. This is the closest one can go to 
// hbook and ID based histogramming
// -------------------------------------------------------------------------
template <class T>
Bool_t AnaBase::fillHist1D(const string& hname, T value, double w) 
{
  TObject *obj = gDirectory->GetList()->FindObject(hname.c_str()); 
  TH1 *h = 0;

  if (obj->InheritsFrom("TH1D"))
    h = dynamic_cast<TH1D*>(obj);
  else if (obj->InheritsFrom("TH1C"))
    h = dynamic_cast<TH1C*>(obj);
  else if (obj->InheritsFrom("TH1K"))
    h = dynamic_cast<TH1K*>(obj);
  else if (obj->InheritsFrom("TH1S"))
    h = dynamic_cast<TH1S*>(obj);
  else if (obj->InheritsFrom("TH1I"))
    h = dynamic_cast<TH1I*>(obj);
  else
    h = dynamic_cast<TH1F*>(obj);

  if (!h) {
    cerr << "**** 1D Histogram <<" << hname << ">> not found" << endl;
    return kFALSE;
  }
  h->Fill(value, w);
  return kTRUE;
}
// ---------------------------------------------
// Convenience routine for filling 2D histograms
// ---------------------------------------------
template <class T1, class T2>
Bool_t AnaBase::fillHist2D(const string& hname, T1 xvalue, T2 yvalue, double w) 
{
  TObject *obj = gDirectory->GetList()->FindObject(hname.c_str()); 
  TH2 *h = 0;

  if (obj->InheritsFrom("TH2D"))
    h = dynamic_cast<TH2D*>(obj);
  else if (obj->InheritsFrom("TH2C"))
    h = dynamic_cast<TH2C*>(obj);
  else if (obj->InheritsFrom("TH2S"))
    h = dynamic_cast<TH2S*>(obj);
  else if (obj->InheritsFrom("TH2I"))
    h = dynamic_cast<TH2I*>(obj);
  else
    h = dynamic_cast<TH2F*>(obj);

  if (!h) {
    cerr << "**** 2D Histogram <<" << hname << ">> not found" << endl;
    return kFALSE;
  }
  h->Fill(xvalue, yvalue, w);
  return kTRUE;
}
// --------------------------------------------------
// Convenience routine for filling profile histograms
// --------------------------------------------------
bool AnaBase::fillProfile(const string& hname, float xvalue, float yvalue, double w) 
{
  TProfile *h = dynamic_cast<TProfile*>(gDirectory->GetList()->FindObject(hname.c_str()));
  if (!h) {
    cerr << "**** Profile Histogram <<" << hname << ">> not found" << endl;
    return false;
  }
  h->Fill(xvalue, yvalue, w);
  return true;
}
// ------------------------------------------------------
// Open the output file with a global filehandle, C++ way
// ------------------------------------------------------
bool AnaBase::openFiles() 
{
  _fLog.open(_logFile.c_str(), ios::out);
  if (!_fLog) {
    cerr << "File: " << _logFile << " could not be opened!" << endl;
    return false;
  }
  _fLog << setiosflags(ios::fixed);
  return true;
}
// ------------------------
// Close the output file
// ------------------------
void AnaBase::closeFiles() 
{
  if (_fLog) {
    _fLog << resetiosflags(ios::fixed); 
    _fLog.close();
  }
}

void AnaBase::setAddresses() 
{
  map<string, TClonesArray**> mapa;
  mapa.insert ( pair<string, TClonesArray**>("Event", &eventA));
  mapa.insert ( pair<string, TClonesArray**>("Trigger", &triggerA));
  mapa.insert ( pair<string, TClonesArray**>("Vertex", &vertexA));
  mapa.insert ( pair<string, TClonesArray**>("Tau", &tauA));
  mapa.insert ( pair<string, TClonesArray**>("Electron", &electronA));
  mapa.insert ( pair<string, TClonesArray**>("Muon", &muonA));
  mapa.insert ( pair<string, TClonesArray**>("Jet", &jetA));
  mapa.insert ( pair<string, TClonesArray**>("MET", &metA));
  mapa.insert ( pair<string, TClonesArray**>("GenParticle", &genParticleA));
  mapa.insert ( pair<string, TClonesArray**>("GenJet", &genJetA));

  map<string, int*> mapb;
  mapb.insert ( pair<string, int*>("nVertex", &n_vertex));
  mapb.insert ( pair<string, int*>("nTau", &n_tau));
  mapb.insert ( pair<string, int*>("nElectron", &n_electron));
  mapb.insert ( pair<string, int*>("nMuon", &n_muon));
  mapb.insert ( pair<string, int*>("nJet", &n_jet));
  mapb.insert ( pair<string, int*>("nMET", &n_met));
  mapb.insert ( pair<string, int*>("nGenParticle", &n_genparticle));
  mapb.insert ( pair<string, int*>("nGenJet", &n_genjet));

  for (map<string, TClonesArray**>::const_iterator it  = mapa.begin(); 
                                                   it != mapa.end(); ++it)  
  {
    const string& bname = it->first;
    TClonesArray** a = it->second;

    TBranch* branch = _chain->GetBranch(bname.c_str());  // Get branch pointer
    assert(branch);
    _chain->SetBranchAddress(bname.c_str(), a);         // Set branch
  }  

  for (map<string, int*>::const_iterator it  = mapb.begin(); 
                                         it != mapb.end(); ++it)  
  {
    const string& bname = it->first;
    int* a = it->second;

    TBranch* branch = _chain->GetBranch(bname.c_str());  // Get branch pointer
    assert(branch);
    _chain->SetBranchAddress(bname.c_str(), a);         // Set branch
  }  
}
int AnaBase::getEntry(int lflag) const
{
  static const string bnames[] = 
  {
    "Event",
    "Trigger",
    "Vertex", 
    "nVertex",
    "Tau", 
    "nTau",
    "Electron", 
    "nElectron",
    "Muon", 
    "nMuon",
    "Jet", 
    "nJet",
    "MET", 
    "nMET",
    "GenParticle",
    "nGenParticle",
    "GenJet",
    "nGenJet"
  };
  int nbytes = 0;
  for (unsigned int i = 0; i < NEL(bnames); ++i) {
    TBranch* branch = _chain->GetBranch(bnames[i].c_str());
    assert(branch);
    nbytes += branch->GetEntry(lflag);
  }

  return nbytes;
}
// not used yet
void AnaBase::enableBranches() 
{
  static const string bnames[] = 
  {
    "Event",
    "Trigger",
    "Vertex", 
    "nVertex",
    "Tau", 
    "nTau",
    "Electron", 
    "nElectron",
    "Muon", 
    "nMuon",
    "Jet", 
    "nJet",
    "MET", 
    "nMET",
    "GenParticle",
    "nGenParticle",
    "GenJet",
    "nGenJet"
  };
  _chain->SetBranchStatus("*", kFALSE); // Disable all branches
  for (unsigned int i = 0; i < NEL(bnames); ++i) {
    _chain->SetBranchStatus(bnames[i].c_str(), kTRUE);
  }
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
bool AnaBase::readJob(const string& jobFile, int& nFiles)
{
  // Define parameters and set default values
  vector<string> fileList;
  string outputFile("pippo.out");

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
  hmap.insert ( pair<string, map<string, double>* >("vtxCutList", &_vtxCutMap));
  hmap.insert ( pair<string, map<string, double>* >("electronCutList", &_electronCutMap));
  hmap.insert ( pair<string, map<string, double>* >("muonCutList", &_muonCutMap));
  hmap.insert ( pair<string, map<string, double>* >("tauCutList", &_tauCutMap));
  hmap.insert ( pair<string, map<string, double>* >("bjetCutList", &_bjetCutMap));

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
    if (key == "logFile")
      _logFile  = tokens[1];
    else if (key == "logOption") 
      _logOption = strtol(tokens[1].c_str(), NULL, 2);
    else if (key == "maxEvent") 
      _maxEvt = atoi(tokens[1].c_str());
    else if (key == "histFile") 
      _histFile = tokens[1];
    else if (key == "inputFile")  {
      // Treat the 'inputFile' key specially
      // In case the input ntuples files follow a pattern like file_[n].root
      // we can specify the input file in a much simpler way, like
      // inputFile start end_index filename omitting '[n][.root]'
      // Treat the 'inputFile' key specially
      if (size >= 3) {
	char name[256];
	int start = (size == 3) ? 1 : atoi(tokens[1].c_str());
	int   end = (size == 3) ? atoi(tokens[1].c_str()) : atoi(tokens[2].c_str());
	for (int i = start; i <= end; i++) {
	  // assumes the file name ends with .root
	  // e.g topnt_[version]_[dataset]_[n]-[stream].root
	  // We need to specify [stream]
          sprintf(name, "%s%d.root", tokens[3].c_str(), i); // assumes the file name ends with .root
	  fileList.push_back(name);
	}
      }
      else if (size == 2) {
	fileList.push_back(tokens[1]);
      }
      else {
	cerr << "**** Error in input Root filenames!" << endl;
	exit(-1);
      }
    }
    else {
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
    }
    tokens.clear();
  }
  // Close the file
  fin.close();

  nFiles = fileList.size();
  if (!nFiles) {
    cerr << "Input Root file list is empty! exiting ..." << endl;
    return false;
  }
  // Build the chain of root files
  for (unsigned int i = 0; i < fileList.size(); i++)
    setInputFile(fileList[i].c_str());

  printJob(fileList);
  return true;
}
void AnaBase::printJob(vector<string>& fileList, ostream& os)
{
  map<string, map<string, double> > hmap;
  hmap.insert ( pair<string, map<string, double> >("vtxCutList", _vtxCutMap));
  hmap.insert ( pair<string, map<string, double> >("electronCutList", _electronCutMap));
  hmap.insert ( pair<string, map<string, double> >("muonCutList", _muonCutMap));
  hmap.insert ( pair<string, map<string, double> >("tauCutList", _tauCutMap));
  hmap.insert ( pair<string, map<string, double> >("bjetCutList", _bjetCutMap));

  os << "       logFile = " <<  _logFile << endl 
     << "      histFile = " <<  _histFile << endl
     << "     logOption = " << _logOption << endl
     << "     maxEvent = " << _maxEvt << endl;

  os << "inputFiles: Total # = " << fileList.size() << endl; 
  for (vector<string>::const_iterator it  = fileList.begin(); 
       it != fileList.end(); ++it) 
    os << *it << endl; 

  for (map<string, map<string, double> >::const_iterator it  = hmap.begin(); 
                                                         it != hmap.end(); ++it)  
  {
    os << "=>>> " << it->first << endl; 
    map<string, double> m = it->second;
    os << setprecision(2);
    for (map<string,double>::const_iterator jt  = m.begin(); 
                                            jt != m.end(); ++jt)  
      os << jt->first << ": " << setw(7) << jt->second << endl;
    os << endl; 
  }
}
