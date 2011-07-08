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

using namespace std;
// -----------
// Constructor
// -----------
AnaBase::AnaBase(const string& filename)
  : _chain(new TChain("treeCreator/vhtree")),
    _histf(0),
          eventA(new TClonesArray("Event")),
       //genEventA(new TClonesArray("GenEvent")),   
       electronA(new TClonesArray("Electron")),
    genParticleA(new TClonesArray("GenParticle")),  
         getJetA(new TClonesArray("GenJet")),   
         genMETA(new TClonesArray("GenMET")),
            metA(new TClonesArray("MET")),    
            tauA(new TClonesArray("Tau")),    
           muonA(new TClonesArray("Muon")),    
            jetA(new TClonesArray("Jet")),    
         vertexA(new TClonesArray("Vertex")),    
        triggerA(new TClonesArray("Trigger"))    
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
//  delete genEventA;
  delete electronA;
  delete genParticleA;
  delete getJetA;
  delete genMETA;
  delete metA;
  delete tauA;
  delete muonA;
  delete jetA;
  delete vertexA;
  delete triggerA;
}
// ------------------------
// Clear the clones arrays
// ------------------------
void AnaBase::clearEvent() 
{
  if (eventA) eventA->Clear();
//  if (genEventA) genEventA->Clear();
  if (electronA) electronA->Clear();
  if (genParticleA) genParticleA->Clear();
  if (getJetA) getJetA->Clear();
  if (genMETA) genMETA->Clear();
  if (metA) metA->Clear();
  if (tauA) tauA->Clear();
  if (muonA) muonA->Clear();
  if (jetA) jetA->Clear();
  if (vertexA) vertexA->Clear();
  if (triggerA) triggerA->Clear();

  n_electron = 0;
  n_jet  = 0;
  n_muon  = 0;
  n_met  = 0;
  n_tau = 0;
  n_vertex = 0;  
  n_genjet = 0;
  n_genmet = 0;
  n_genparticle = 0;
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool AnaBase::beginJob() 
{
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
}
// -------------------
// The main event loop
// -------------------
void AnaBase::eventLoop() 
{
  const int nPrint = 1000;

  // Initialize analysis
  if (!beginJob()) return;

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  for (int ev = 0; ev < nEvents; ev++) {
    bool select = false;

    clearEvent();
    int lflag = _chain->LoadTree(ev); 
    int nentries = getEntry(lflag);    // returns bytes read
    
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
    if (ev%nPrint == 0) cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
                             << " ==> " << _chain->GetCurrentFile()->GetName() 
                             << " <<< Run# " << run 
                             << " Event# " << setw(8) << event << " >>> " 
                             << " Events proc. " << setw(8) << ev << endl;

    // Let's look at the Tau collection
    _fLog << setprecision(3);
    if (n_tau) {
      _fLog << "=>> Taus: " << n_tau << endl;
      _fLog << "indx     Eta     Phi      Pt  Energy"
  	    << "   lchPt   lnpPt    lpPt" 
	    << endl; 
      for (int indx = 0; indx < n_tau; ++indx) {
        const Tau* tau = dynamic_cast<Tau*>(tauA->At(indx));
        if (!tau) continue;
        if (abs(tau->eta) >= _tauCutMap["eta"] ||
              tau->decayModeFinding <= 0.5 ||            
              tau->looseIsolation <= 0.5 ||            
              tau->againstMuonTight <= 0.5 ||            
              tau->againstElectronLoose <= 0.5
	    ) continue;
        _fLog << setprecision(2)
              << setw(4) << indx 
              << setw(8) << tau->eta
              << setw(8) << tau->phi
              << setw(8) << tau->pt
              << setw(8) << tau->energy
              << setw(8) << tau->leadChargedParticlePt
              << setw(8) << tau->leadNeutralParticlePt
              << setw(8) << tau->leadParticlePt
              << endl;
      }
    }
    // Let's look at the Muon collection
    if (n_muon) {
      _fLog << "=>> Muons: " << n_muon << endl;
      _fLog << "indx     Eta     Phi      Pt       P"
  	    << "      D0   D0Err      Dz   DzErr" 
	    << endl; 
      for (int indx = 0; indx < n_muon; ++indx) {
        const Muon* muon = dynamic_cast<Muon*>(muonA->At(indx));
        if (!muon) continue;
        // Now apply cuts
        if (!muon->isTrackerMuon                  ||
            fabs(muon->eta) >= _muonCutMap["eta"] || 
            muon->relIso >= _muonCutMap["relIso"] ||
            (muon->pixHits + muon->trkHits) <= _muonCutMap["ptHits"] ||
            muon->globalChi2 >= _muonCutMap["globalChi2"] ||
            fabs(muon->trkD0) >= _muonCutMap["trkD0"]
	   ) continue;
        _fLog << setw(4) << indx 
              << setw(8) << muon->eta
              << setw(8) << muon->phi
              << setw(8) << muon->pt
              << setw(8) << muon->p
              << setw(8) << muon->trkD0
              << setw(8) << muon->trkD0Error
              << setw(8) << muon->trkDz
              << setw(8) << muon->trkDzError
              << endl;
	fillHist1D("muonPt", muon->pt, 1.0);
      }
    }
  }  
  // Analysis is over
  endJob();
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void AnaBase::endJob() 
{
  cout << resetiosflags(ios::fixed);

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
  unsigned int found = fname.find("root:");
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
  TBranch* branch = _chain->GetBranch("Event");     // Get branch pointer
  assert(branch);
  _chain->SetBranchAddress("Event", &eventA);       // Set b

  branch = _chain->GetBranch("Trigger");
  assert(branch);
  _chain->SetBranchAddress("Trigger", &triggerA);

  branch = _chain->GetBranch("Vertex");
  assert(branch);
  _chain->SetBranchAddress("Vertex", &vertexA);

  branch = _chain->GetBranch("nTau");
  assert(branch);
  _chain->SetBranchAddress("nTau", &n_tau);

  branch = _chain->GetBranch("Tau");
  assert(branch);
  _chain->SetBranchAddress("Tau", &tauA);

  branch = _chain->GetBranch("nElectron");
  assert(branch);
  _chain->SetBranchAddress("nElectron", &n_electron);

  branch = _chain->GetBranch("Electron");
  assert(branch);
  _chain->SetBranchAddress("Electron", &electronA);

  branch = _chain->GetBranch("nMuon");
  assert(branch);
  _chain->SetBranchAddress("nMuon", &n_muon);

  branch = _chain->GetBranch("Muon");
  assert(branch);
  _chain->SetBranchAddress("Muon", &muonA);

  branch = _chain->GetBranch("nJet");
  assert(branch);
  _chain->SetBranchAddress("nJet", &n_jet);

  branch = _chain->GetBranch("Jet");
  assert(branch);
  _chain->SetBranchAddress("Jet", &jetA);
}
int AnaBase::getEntry(int lflag) const
{
  TBranch* branch = _chain->GetBranch("Event");
  assert(branch);
  int nbytes = branch->GetEntry(lflag);

  branch = _chain->GetBranch("Trigger");
  assert(branch);
  nbytes += branch->GetEntry(lflag);

  branch = _chain->GetBranch("Vertex");
  assert(branch);
  nbytes += branch->GetEntry(lflag);

  branch = _chain->GetBranch("nVertex");
  assert(branch);
  nbytes += branch->GetEntry(lflag);

  branch = _chain->GetBranch("Electron");
  assert(branch);
  nbytes += branch->GetEntry(lflag);

  branch = _chain->GetBranch("nElectron");
  assert(branch);
  nbytes += branch->GetEntry(lflag);

  branch = _chain->GetBranch("Tau");
  assert(branch);
  nbytes += branch->GetEntry(lflag);

  branch = _chain->GetBranch("nTau");
  assert(branch);
  nbytes += branch->GetEntry(lflag);

  branch = _chain->GetBranch("Muon");
  assert(branch);
  nbytes += branch->GetEntry(lflag);

  branch = _chain->GetBranch("nMuon");
  assert(branch);
  nbytes += branch->GetEntry(lflag);

  return nbytes;
}
void AnaBase::enableBranches() 
{
  _chain->SetBranchStatus("*", kFALSE); // Disable all branches

  // Now enable one by one
  _chain->SetBranchStatus("Event", kTRUE);
//  _chain->SetBranchStatus("GenEvent", kTRUE);
  _chain->SetBranchStatus("Trigger", kTRUE);
  _chain->SetBranchStatus("Vertex", kTRUE);
  _chain->SetBranchStatus("nVertex", kTRUE);
  _chain->SetBranchStatus("Electron", kTRUE);
  _chain->SetBranchStatus("nElectron", kTRUE);
  _chain->SetBranchStatus("Muon", kTRUE);
  _chain->SetBranchStatus("nMuon", kTRUE);
  _chain->SetBranchStatus("Tau", kTRUE);
  _chain->SetBranchStatus("nTau", kTRUE);
  _chain->SetBranchStatus("MET", kTRUE);
  _chain->SetBranchStatus("nMET", kTRUE);
  _chain->SetBranchStatus("Jet", kTRUE);
  _chain->SetBranchStatus("nJet", kTRUE);
  _chain->SetBranchStatus("GenParticle", kTRUE);
  _chain->SetBranchStatus("nGenParticle", kTRUE);
  _chain->SetBranchStatus("GenJet", kTRUE);
  _chain->SetBranchStatus("nGenJet", kTRUE);
  _chain->SetBranchStatus("GenParticle", kTRUE);
  _chain->SetBranchStatus("nGenParticle", kTRUE);
  _chain->SetBranchStatus("GenParticle", kTRUE);
  _chain->SetBranchStatus("nGenParticle", kTRUE);
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
    else if (key == "muonCutList") {
      // Muon Cuts
      _muonCutMap.clear();
      for (int i = 1; i < size; ++i) {
        vector<string> cutstr;
        // Split the line into words
        AnaUtil::tokenize(tokens[i], cutstr, "=");
        _muonCutMap.insert ( pair<string,double>(cutstr[0], atof(cutstr[1].c_str())));
      }
    }
    else if (key == "tauCutList") {
      // Tau Cuts
      _tauCutMap.clear();
      for (int i = 1; i < size; ++i) {
        vector<string> cutstr;
        // Split the line into words
        AnaUtil::tokenize(tokens[i], cutstr, "=");
        _tauCutMap.insert ( pair<string,double>(cutstr[0], atof(cutstr[1].c_str())));
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
  os << "       logFile = " <<  _logFile << endl 
     << "      histFile = " <<  _histFile << endl;

  os << "inputFiles: Total # = " << fileList.size() << endl; 
  for (vector<string>::const_iterator it  = fileList.begin(); 
                                      it != fileList.end(); ++it) 
    os << *it << endl; 

  os << "=>>> muonCutList: " << endl; 
  os << setprecision(2);
  for (map<string,double> ::const_iterator it  = _muonCutMap.begin(); 
                                           it != _muonCutMap.end(); ++it)  
    os << (*it).first << " => " << setw(7) << (*it).second << endl;
  os << endl; 

  os << "=>>> tauCutList: " << endl; 
  for (map<string,double> ::const_iterator it  = _tauCutMap.begin(); 
                                           it != _tauCutMap.end(); ++it)  
    os << (*it).first << " => " << setw(7) << (*it).second << endl;
  os << endl; 
}
