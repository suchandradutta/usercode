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
#include <cmath>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TProfile.h"

#include "AnaBase.h"
#include "AnaUtil.h"

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
AnaBase::AnaBase()
  : _chain(new TChain("treeCreator/vhtree")),
    _histf(0),
          eventA(new TClonesArray("vhtm::Event")),
         vertexA(new TClonesArray("vhtm::Vertex")),    
            tauA(new TClonesArray("vhtm::Tau")),    
       electronA(new TClonesArray("vhtm::Electron")),
           muonA(new TClonesArray("vhtm::Muon")),    
            jetA(new TClonesArray("vhtm::Jet")),    
            metA(new TClonesArray("vhtm::MET")),    
    genParticleA(new TClonesArray("vhtm::GenParticle")),  
         genJetA(new TClonesArray("vhtm::GenJet")),   
         genMetA(new TClonesArray("vhtm::GenMET")),
     triggerobjA(new TClonesArray("vhtm::TriggerObject")),
          trackA(new TClonesArray("vhtm::Track")),
       //genEventA(new TClonesArray("vhtm::GenEvent")),   
          _l1physbits(new vector<int>()),
          _l1techbits(new vector<int>()),
            _hltpaths(new vector<string>()),
          _hltresults(new vector<int>()),
        _hltprescales(new vector<int>()),
        _isMC(false),
        _readTrk(false),
        _readTrigObject(true),
        _logOption(0),
        _studyTrigger(true),
        _histFile("default.root"),
        _puHistFile("./reweightFunctionFall11.root"),
        _logFile("default.out"),
        _evFile("default_event.out"),
        _maxEvt(0)
{
  cout << setiosflags(ios::fixed); 
  cout << "=== Start of Analysis === " << endl;
  _fileList.clear();
  _brList.clear();
  _puWtList.clear();
  _trigPathList.clear();
}
// ----------
// Destructor
// ----------
AnaBase::~AnaBase() 
{
  clearEvent();

  delete eventA;
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
  delete triggerobjA;
  delete trackA;
}
// ------------------------
// Clear the clones arrays
// ------------------------
void AnaBase::clearEvent() 
{
  if (eventA) eventA->Clear();
//  if (genEventA) genEventA->Clear();
  if (vertexA) vertexA->Clear();
  if (tauA) tauA->Clear();
  if (electronA) electronA->Clear();
  if (muonA) muonA->Clear();
  if (jetA) jetA->Clear();
  if (metA) metA->Clear();
  if (genParticleA) genParticleA->Clear();
  if (genJetA) genJetA->Clear();
  if (genMetA) genMetA->Clear();
  if (triggerobjA) triggerobjA->Clear();
  if (trackA) trackA->Clear();

  if (_l1physbits)   _l1physbits->clear();
  if (_l1techbits)   _l1techbits->clear();
  if (_hltpaths)     _hltpaths->clear();
  if (_hltresults)   _hltresults->clear();
  if (_hltprescales) _hltprescales->clear();

  n_vertex      = 0;  
  n_tau         = 0;
  n_electron    = 0;
  n_muon        = 0;
  n_jet         = 0;
  n_met         = 0;
  n_genparticle = 0;
  n_genjet      = 0;
  n_genmet      = 0;
  n_triggerobj  = 0;
  n_track       = 0;
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool AnaBase::beginJob() 
{ 
  if (_isMC) readPileUpHist();

  // Open the output ROOT file
  _histf = TFile::Open(_histFile.c_str(), "RECREATE");

  setAddresses();
  //enableBranches();
  nEvents = static_cast<int>(_chain->GetEntries()); 
  if (nEvents <= 0) {
    cerr << "******* nEvents = " << nEvents << ", returning!" << endl;
    return false;
  }
  if (_maxEvt > 0) nEvents = std::min(nEvents,_maxEvt);
  cout << " ===== # of events to analyse, nEvents = " << nEvents << endl;

  openFiles();

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
int AnaBase::setInputFile(const string& fname) 
{
  size_t found = fname.find("root:");
  if (found == string::npos && gSystem->AccessPathName(fname.c_str())) {
     cerr << ">>> Warning: File <<" << fname << ">> was not found!!" << endl;
     return static_cast<int>(_chain->GetEntries()); 
  }
  _chain->AddFile(fname.c_str(), -1);
  return static_cast<int>(_chain->GetEntries()); 
}
// ---------------------------------------
// Get total number of events in the chain
// --------------------------------------
int AnaBase::getEntries() const 
{
  return static_cast<int>(_chain->GetEntries());
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

  _evLog.open(_evFile.c_str(), ios::out);
  if (!_evLog) {
    cerr << "File: " << _evFile << " could not be opened!" << endl;
    return false;
  }
  _evLog << setiosflags(ios::fixed);
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
  if (_evLog) {
    _evLog << resetiosflags(ios::fixed); 
    _evLog.close();
  }
}
void AnaBase::setAddresses() 
{
  map<string, TClonesArray**> mapa;
  mapa.insert(pair<string, TClonesArray**>("Event", &eventA));
  mapa.insert(pair<string, TClonesArray**>("Vertex", &vertexA));
  mapa.insert(pair<string, TClonesArray**>("Tau", &tauA));
  mapa.insert(pair<string, TClonesArray**>("Electron", &electronA));
  mapa.insert(pair<string, TClonesArray**>("Muon", &muonA));
  mapa.insert(pair<string, TClonesArray**>("Jet", &jetA));
  mapa.insert(pair<string, TClonesArray**>("MET", &metA));
  if (_readTrigObject) mapa.insert(pair<string, TClonesArray**>("TriggerObject", &triggerobjA));
  if (_readTrk) mapa.insert(pair<string, TClonesArray**>("Track", &trackA));
  if (_isMC) {
    mapa.insert(pair<string, TClonesArray**>("GenParticle", &genParticleA));
    mapa.insert(pair<string, TClonesArray**>("GenJet", &genJetA));
  }

  map<string, int*> mapb;
  mapb.insert(pair<string, int*>("nVertex", &n_vertex));
  mapb.insert(pair<string, int*>("nTau", &n_tau));
  mapb.insert(pair<string, int*>("nElectron", &n_electron));
  mapb.insert(pair<string, int*>("nMuon", &n_muon));
  mapb.insert(pair<string, int*>("nJet", &n_jet));
  mapb.insert(pair<string, int*>("nMET", &n_met));
  if (_readTrigObject) mapb.insert(pair<string, int*>("nTriggerObject", &n_triggerobj));
  if (_readTrk) mapb.insert(pair<string, int*>("nTrack", &n_track));
  if (_isMC) {
    mapb.insert(pair<string, int*>("nGenParticle", &n_genparticle));
    mapb.insert(pair<string, int*>("nGenJet", &n_genjet));
  }

  for (map<string, TClonesArray**>::const_iterator it  = mapa.begin(); 
                                                   it != mapa.end(); ++it)  
  {
    const string& bname = it->first;
    TClonesArray** a = it->second;
    TBranch* branch = _chain->GetBranch(bname.c_str());  // Get branch pointer
    if (!branch) {
       cout << ">>> SetBranchAddress: <" << bname << "> not found!" << endl;
       continue;
    }
    cout << ">>> SetBranchAddress: <" << bname << ">"  << endl;
    _chain->SetBranchAddress(bname.c_str(), a);         // Set branch
    _brList.push_back(bname);
  }  

  for (map<string, int*>::const_iterator it  = mapb.begin(); 
                                         it != mapb.end(); ++it)  
  {
    const string& bname = it->first;
    int* a = it->second;

    TBranch* branch = _chain->GetBranch(bname.c_str());  // Get branch pointer
    if (!branch) {
       cout << ">>> SetBranchAddress: <" << bname << "> not found!" << endl;
       continue;
    }
    cout << ">>> SetBranchAddress: <" << bname << ">"  << endl;
    _chain->SetBranchAddress(bname.c_str(), a);         // Set branch
    _brList.push_back(bname);
  }  

  // Now the trigger variables
  map<string, vector<int>**> mapc;
  mapc.insert(pair<string, vector<int>**>("l1physbits", &_l1physbits));
  mapc.insert(pair<string, vector<int>**>("l1techbits", &_l1techbits));
  mapc.insert(pair<string, vector<int>**>("hltresults", &_hltresults));
  mapc.insert(pair<string, vector<int>**>("hltprescales", &_hltprescales));

  for (map<string, vector<int>**>::const_iterator it  = mapc.begin(); 
                                                  it != mapc.end(); ++it)  
  {
    const string& bname = it->first;
    vector<int>** a = it->second;

    TBranch* branch = _chain->GetBranch(bname.c_str());  // Get branch pointer
    if (!branch) {
       cout << ">>> SetBranchAddress: <" << bname << "> not found!" << endl;
       continue;
    }
    cout << ">>> SetBranchAddress: <" << bname << ">"  << endl;
   _chain->SetBranchAddress(bname.c_str(), a);         // Set branch
    _brList.push_back(bname);
  }  

  TBranch* branch = _chain->GetBranch("hltpaths");  // Get branch pointer                                                                                  
  if (!branch) {
     cout << ">>> SetBranchAddress: <hltpaths> not found!" << endl;
     return;
  }
  cout << ">>> SetBranchAddress: <hltpaths>"  << endl;
  _chain->SetBranchAddress("hltpaths", &_hltpaths);
  _brList.push_back("hltpaths");
}
int AnaBase::getEntry(int lflag) const
{
  int nbytes = 0;
  for (vector<string>::const_iterator it  = _brList.begin(); 
                                      it != _brList.end(); ++it) {
    TBranch* branch = _chain->GetBranch((*it).c_str());
    if (!branch) {
       cout << ">>> Branch: " << (*it) << " not found!" << endl;
       continue;
    }
    nbytes += branch->GetEntry(lflag);
  }
  return nbytes;
}
// not used yet
void AnaBase::enableBranches() 
{
  _chain->SetBranchStatus("*", kFALSE); // Disable all branches
  for (vector<string>::const_iterator it  = _brList.begin(); 
                                      it != _brList.end(); ++it) {
    _chain->SetBranchStatus((*it).c_str(), kTRUE);
  }
}
bool AnaBase::readJob(const string& jobFile, int& nFiles)
{
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
    string key = tokens.at(0);
    string value = tokens.at(1);
    if (key == "dataType") 
      _isMC = (value == "mc" || value == "MC") ? true : false;
    else if (key == "readTrk") 
      _readTrk = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "readTrigObject") 
      _readTrigObject = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "studyTrigger") 
      _studyTrigger = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "logFile")
      _logFile = value;
    else if (key == "eventFile")
      _evFile  = value;
    else if (key == "logOption") 
      _logOption = strtol(value.c_str(), NULL, 2);
    else if (key == "maxEvent") 
      _maxEvt = atoi(value.c_str());
    else if (key == "histFile") 
      _histFile = value;
    else if (key == "puHistFile") 
      _puHistFile = value;
    else if (key == "trigPathList") {
      for (size_t i = 1; i < tokens.size(); ++i) {
	_trigPathList.push_back(tokens.at(i));       
      }
    }
    else if (key == "inputFile")  {
      // Treat the 'inputFile' key specially
      // In case the input ntuples files follow a pattern like file_[n].root
      // we can specify the input file in a much simpler way, like
      // inputFile start end_index filename omitting '[n][.root]'
      // Treat the 'inputFile' key specially
      if (size >= 3) {
	char name[256];
	int start = (size == 3) ? 1 : atoi(value.c_str());
	int   end = (size == 3) ? atoi(value.c_str()) : atoi(tokens[2].c_str());
	for (int i = start; i <= end; i++) {
	  // assumes the file name ends with .root
	  // e.g topnt_[version]_[dataset]_[n]-[stream].root
	  // We need to specify [stream]
          sprintf(name, "%s%d.root", tokens[3].c_str(), i); // assumes the file name ends with .root
	  _fileList.push_back(name);
	}
      }
      else if (size == 2) {
	_fileList.push_back(value);
      }
      else {
	cerr << "**** Error in input Root filenames!" << endl;
	exit(-1);
      }
    }
    tokens.clear();
  }
  // Close the file
  fin.close();

  // Build the chain of root files
  for (vector<string>::const_iterator it  = _fileList.begin();
                                      it != _fileList.end(); ++it) {
     string fname = *it;
     cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
     ++nFiles;
     int nevt = setInputFile(fname.c_str());
     if (_maxEvt > 0 && nevt >= _maxEvt) break;
  }

  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  }

  return true;
}
void AnaBase::printJob(ostream& os) const
{
  os << "      datatype = " << ((_isMC) ? "mc" : "data") << endl
     << "       logFile = " << _logFile << endl 
     << "     eventFile = " << _evFile << endl
     << "      histFile = " << _histFile << endl
     << "    puHistFile = " << _puHistFile << endl
     << "  studyTrigger = " << _studyTrigger << endl
     << "     logOption = " << _logOption << endl
     << "      maxEvent = " << _maxEvt << endl;

  // Trigger Path List
  os << ">>> INFO. Trigger Paths used:" << endl;
  for (vector<string>::const_iterator it  = _trigPathList.begin();
                                      it != _trigPathList.end(); ++it)
    os << *it << endl;

  // InputFiles
  if (_chain) {
     TObjArray *fileElements = _chain->GetListOfFiles();
     os << ">>> INFO. nFiles: " << fileElements->GetEntries() 
        << ", Files to analyse:" << endl;
     TIter next(fileElements);
     TChainElement *chEl = 0;
     while (( chEl = dynamic_cast<TChainElement*>(next()) ))
       os << chEl->GetTitle() << endl;
  }
  else {
     os << ">>> INFO. inputFiles: Total # = " << _fileList.size() << endl;
     for (vector<string>::const_iterator it  = _fileList.begin();
                                         it != _fileList.end(); ++it)
       os << *it << endl;
  }
}
// Collect Object information
// Let's look at the event vertex  
void AnaBase::findVtxInfo(vector<Vertex>& list, map<string, double> cutMap, Options& op, ostream& os) const {
  if (n_vertex < 1) return;

  if (op.verbose)
    os << "=>> Vertices: " << n_vertex << endl
       << "indx     ndf     dxy       z   sumPt    chi2   ntrks  ntrkw05     sbit" 
       << endl; 
  for (int indx = 0; indx < n_vertex; ++indx) {
    const Vertex* vtx = dynamic_cast<Vertex*>(vertexA->At(indx));
    if (!vtx) continue;

    double dxy = sqrt(pow(vtx->x, 2) + pow(vtx->y, 2));
    int sbit = 0;
    if (vtx->ndf <= AnaUtil::cutValue(cutMap, "ndf"))    sbit |= (1 << 0);
    if (dxy >= AnaUtil::cutValue(cutMap, "dxy"))         sbit |= (1 << 1);
    if (abs(vtx->z) >= AnaUtil::cutValue(cutMap, "z"))   sbit |= (1 << 2);

    if (op.verbose) {
      os << setprecision(2)
         << setw(4) << indx
         << setw(8) << vtx->ndf
         << setw(8) << dxy
         << setw(8) << vtx->z 
         << setw(8) << vtx->sumPt
         << setw(8) << vtx->chi2
         << setw(8) << vtx->ntracks
         << setw(9) << vtx->ntracksw05;
      AnaUtil::bit_print(sbit, 8, os);
    }
    if (op.usesbit && sbit) continue;
    list.push_back(*vtx);
  }
  //if (list.size() > 1) 
  //  sort(list.begin(), list.end(), VertexComparator());
}
// Let's look at the Muon collection
void AnaBase::findMuonInfo(vector<Muon>& list,  map<string, double> cutMap, Options& op, ostream& os) const {
  if (n_muon < 1) return;

  if (op.verbose)
    os << "=>> Muons: " << n_muon << endl
          << "indx     Eta     Phi      Pt       P"
          << "      D0   D0Err      Dz   DzErr" 
          << "   relIso  pfRelIso pixHits trkHits   gChi2      dB            sbit"
          << endl; 

  for (int indx = 0; indx < n_muon; ++indx) {
    const Muon* muon = dynamic_cast<Muon*>(muonA->At(indx));
    if (!muon) continue;

    int sbit = 0;
    if (!muon->isTrackerMuon)                                                   sbit |= (1 <<  0);
    if (!muon->isGlobalMuonPromptTight)                                         sbit |= (1 <<  1);
    if (!muon->isAllArbitrated)                                                 sbit |= (1 <<  2);
    if (abs(muon->vtxDistZ) >= AnaUtil::cutValue(cutMap, "vtxDistZ"))           sbit |= (1 <<  3);
    if (muon->nChambers < AnaUtil::cutValue(cutMap, "nChambers"))               sbit |= (1 <<  4);
    if (muon->nMatches < AnaUtil::cutValue(cutMap, "nMatches"))                 sbit |= (1 <<  5);
    if (muon->nMatchedStations < AnaUtil::cutValue(cutMap, "nMatchedStations")) sbit |= (1 <<  6); 
    if (abs(muon->eta) >= AnaUtil::cutValue(cutMap, "eta"))                     sbit |= (1 <<  7);
    if (muon->pfRelIso >= AnaUtil::cutValue(cutMap, "relIso"))                  sbit |= (1 <<  8);
    if (muon->trkHits < AnaUtil::cutValue(cutMap,"trkHits"))                    sbit |= (1 <<  9);
    if (muon->pixHits < AnaUtil::cutValue(cutMap,"pixHits"))                    sbit |= (1 << 10);
    if (muon->globalChi2 >= AnaUtil::cutValue(cutMap,"globalChi2"))             sbit |= (1 << 11);
    if (abs(muon->trkD0) >= AnaUtil::cutValue(cutMap,"trkD0"))                  sbit |= (1 << 12);
    if (abs(muon->dB) >= AnaUtil::cutValue(cutMap,"dB"))                        sbit |= (1 << 13);

    if (op.verbose) {
      os << setprecision(2)
         << setw(4) << indx 
         << setw(8) << muon->eta
         << setw(8) << muon->phi
         << setw(8) << muon->pt
         << setw(8) << muon->p
         << setw(8) << muon->trkD0
         << setw(8) << muon->trkD0Error
         << setw(8) << muon->trkDz
         << setw(8) << muon->trkDzError
         << setprecision(3)
         << setw(9) << muon->relIso
         << setw(9) << muon->pfRelIso
         << setprecision(3)
         << setw(8) << muon->pixHits
         << setw(8) << muon->trkHits
         << setw(8) << muon->globalChi2
         << setw(8) << muon->dB;
      AnaUtil::bit_print(sbit, 13, os);
    }
    // Now apply cuts
    if (op.usesbit && sbit) continue;
    list.push_back(*muon);
  }
  if (list.size() > 1) 
    sort(list.begin(), list.end(), PtComparator<Muon>());
}
// Let's look at the Electron collection now
void AnaBase::findElectronInfo(vector<Electron>& list, map<string, double> cutMap, Options& op, ostream& os) const {
  if (n_electron < 1) return;
  if (op.verbose)
    os << "=>> Electrons: " << n_electron << endl
       << "indx     Eta     Phi      Pt  Energy"
       << "   scEta   dB      eleId     sbit" 
       << endl; 

  for (int indx = 0; indx < n_electron; ++indx) {
    const Electron* elec = dynamic_cast<Electron*>(electronA->At(indx));
    if (!elec) continue;

    bool quality_EB_loose = (abs(elec->eta) <= 1.4442 
                         &&  elec->sigmaEtaEta < 0.01 
                         &&  elec->deltaEtaTrkSC < 0.007 
                         &&  elec->deltaPhiTrkSC < 0.8
                         &&  elec->hoe < 0.15);
    bool quality_EE_loose = (abs(elec->eta) >= 1.566 
                         &&  elec->sigmaEtaEta < 0.03
                         &&  elec->deltaEtaTrkSC < 0.01 
                         &&  elec->deltaPhiTrkSC < 0.7
                         &&  elec->hoe < 0.07);
    bool quality_loose = quality_EB_loose || quality_EE_loose;  

    int sbit = 0;
    if (elec->pt <= AnaUtil::cutValue(cutMap, "pt"))                             sbit |= (1 << 0);
    if (abs(elec->eta) >= AnaUtil::cutValue(cutMap, "eta"))                      sbit |= (1 << 1);
    if (!quality_loose)                                                          sbit |= (1 << 2);
    if (elec->missingHits > AnaUtil::cutValue(cutMap, "missingHits"))            sbit |= (1 << 3);
    //if (elec->simpleEleId95cIso <= AnaUtil::cutValue(cutMap, "eleId"))           sbit |= (1 << 2);
    //if (!elec->hasGsfTrack)                                                      sbit |= (1 << 3);
    if (abs(elec->dB) >= AnaUtil::cutValue(cutMap, "dB"))                        sbit |= (1 << 4);
    if ( abs(elec->scEta) >= AnaUtil::cutValue(cutMap, "scEtaLow") 
      && abs(elec->scEta) <= AnaUtil::cutValue(cutMap, "scEtaUp"))               sbit |= (1 << 5);

    if (op.verbose) {
      os << setprecision(2)
         << setw(4) << indx 
         << setw(8) << elec->eta
         << setw(8) << elec->phi
         << setw(8) << elec->pt
         << setw(8) << elec->energy
         << setw(8) << elec->scEta
         << setw(8) << elec->dB
         << setw(8) << elec->simpleEleId95cIso;
      AnaUtil::bit_print(sbit, 8, os);
    }
    // Now apply cuts
    if (op.usesbit && sbit) continue;
    list.push_back(*elec);
  }
  if (list.size() > 1) 
    sort(list.begin(), list.end(), PtComparator<Electron>());
}
// Let's look at the Tau collection
void AnaBase::findTauInfo(vector<Tau>& list, map<string, double> cutMap, double vz, Options& op, ostream& os) const {
  if (n_tau < 1) return;
  if (op.verbose)
    os << "=>> Taus: " << n_tau << endl
       << "indx     Eta     Phi      Pt  Energy"
       << "   lchPt   lnpPt    lpPt DMF  LI aMT aEL     sbit" 
       << endl; 
  for (int indx = 0; indx < n_tau; ++indx) {
    const Tau* tau = dynamic_cast<Tau*>(tauA->At(indx));
    if (!tau) continue;

    int sbit = 0;
    if (abs(tau->eta) >= AnaUtil::cutValue(cutMap, "eta"))     sbit |= (1 << 0); 
    if (vz != -999 && abs(tau->zvertex - vz) >= AnaUtil::cutValue(cutMap, "dz")) sbit |= (1 << 1);
    if (tau->decayModeFinding     <= 0.5)                      sbit |= (1 << 2); 
    if (tau->byLooseCombinedIsolationDeltaBetaCorr  <= 0.5)    sbit |= (1 << 3); 
    if (tau->againstMuonTight     <= 0.5)                      sbit |= (1 << 4); 
    if (tau->againstElectronLoose <= 0.5)                      sbit |= (1 << 5); 

    if (op.verbose) {
      os << setprecision(2)
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
         << setw(4) << tau->byLooseCombinedIsolationDeltaBetaCorr
         << setw(4) << tau->againstMuonTight
         << setw(4) << tau->againstElectronLoose;
      AnaUtil::bit_print(sbit, 8, os);
    }
    if (op.usesbit && sbit) continue;
    list.push_back(*tau);
  }
  if (list.size() > 1) 
    sort(list.begin(), list.end(), PtComparator<Tau>());   
}
// Let's look at the Jet collection
void AnaBase::findJetInfo(vector<Jet>& list, map<string, double> cutMap, Options& op, ostream& os) const {
  if (n_jet < 1) return;
  if (op.verbose)
    os << "=>> Jets: " << n_jet << endl
       << "indx     Eta     Phi      Pt  Energy"
       << "    TCHE    TCHP   JPBTag  JBPBTag     sbit"
       << endl; 

  for (int indx = 0; indx < n_jet; ++indx) {
    const Jet* jt = dynamic_cast<Jet*>(jetA->At(indx));
    if (!jt) continue;

    int sbit = 0;
    if (abs(jt->eta) >= AnaUtil::cutValue(cutMap, "eta"))                        sbit |= (1 << 0);
    if (jt->pt <= AnaUtil::cutValue(cutMap, "pt"))                               sbit |= (1 << 1);
    if (jt->trackCountingHighEffBTag <= AnaUtil::cutValue(cutMap, "trackCount")) sbit |= (1 << 2);

    if (op.verbose) {
      os << setprecision(2)
         << setw(4) << indx 
         << setw(8) << jt->eta
         << setw(8) << jt->phi
         << setw(8) << jt->pt
         << setw(8) << jt->energy
         << setprecision(1)
         << setw(8) << jt->trackCountingHighEffBTag  
         << setw(8) << jt->trackCountingHighPurBTag
         << setw(8) << jt->jetProbabilityBTag
         << setw(8) << jt->jetBProbabilityBTag;
      AnaUtil::bit_print(sbit, 8, os);
    }
    // Now apply cuts
    if (op.usesbit && sbit) continue;
    list.push_back(*jt);
  }
  if (list.size() > 1) 
    sort(list.begin(), list.end(), PtComparator<Jet>());
}
void AnaBase::dumpEvent(const vector<map<string, double> >& maps, ostream& os) const {
  // Dump original content present in the tree
  // Event
  const Event* ev = dynamic_cast<Event*>(eventA->At(0));
  assert(ev);
  os << "Event " << ev->event 
     << " Lumis " << ev->lumis 
     << " Run " << ev->run 
     << endl;

  // Options common for all the objects
  Options options;
  options.verbose = true;
  options.usesbit = false;

  vector<Vertex> vList;
  findVtxInfo(vList, maps.at(0), options, os);
  double vz = (vList.size() > 0) ? vList.at(0).z : -999;

  vector<Tau> tList;
  findTauInfo(tList, maps.at(1), vz, options, os);

  vector<Muon> mList;
  findMuonInfo(mList, maps.at(2), options, os);

  vector<Electron> eList;
  findElectronInfo(eList, maps.at(3), options, os);

  vector<Jet> bList;
  findJetInfo(bList, maps.at(4), options, os);
}
void AnaBase::dumpGenInfo(ostream& os) const {
  if (!n_genparticle) return;
  os << setprecision(2);
  os << "indx    status    pdgId     eta      phi      pt     energy             mID                             dID"
        << endl;
  for (int indx = 0; indx < n_genparticle; ++indx) {
    const GenParticle* gp = dynamic_cast<GenParticle*>(genParticleA->At(indx));
    if (!gp) continue;

    std::ostringstream mID;
    vector<int> m = gp->motherIndices;
    for (size_t i = 0; i < m.size(); ++i) {
      int mi = m.at(i);
      if (mi >= n_genparticle) continue;
      const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(mi));
      if (!mgp) continue;
      mID << " " << mgp->pdgId; 
    }
    string ms = mID.str();
    if (!ms.length()) ms = " -";
    
    std::ostringstream dID;
    vector<int> d = gp->daughtIndices;
    for (size_t i = 0; i < d.size(); ++i) {
    	  int di = d.at(i);
      if (di >= n_genparticle) continue;
      const GenParticle* dgp = dynamic_cast<GenParticle*>(genParticleA->At(di));
      if (!dgp) continue;
      double energy = dgp->energy;
      int pdgid = dgp->pdgId;
      if (abs(pdgid) == 21 && energy <= 1) continue;
      dID << " " << dgp->pdgId; 
    }
    string ds = dID.str();
    if (!ds.length()) ds = " -";

    os << setw(4)  << indx
          << setw(8)  << gp->status
          << setw(10) << gp->pdgId
          << setw(10) << gp->eta
          << setw(9)  << gp->phi
          << setw(9)  << gp->pt
          << setw(9)  << gp->energy
          << setw(16) << ms 
          << ds
          << endl;
#if 0 
    vector<int> m = gp->motherIndices;
    os << "# of mother: "  << m.size() << endl;
    for (size_t i = 0; i < m.size(); ++i) {
      int mi = m.at(i);
      os << "index: " << mi << endl;
      if (mi >= n_genparticle) continue;
      const GenParticle* mgp = dynamic_cast<GenParticle*>(genParticleA->At(mi));
      os << "\tpdgId: " << mgp->pdgId << ", status:" << mgp->status << endl; 
    }

    vector<int> d = gp->daughtIndices;
    os << "# of daughter: " << d.size() << endl;
    for (size_t i = 0; i < d.size(); ++i) {
    	  int di = d.at(i);
      os << "index: " << di << endl;
      if (di >= n_genparticle) continue;
      const GenParticle* dgp = dynamic_cast<GenParticle*>(genParticleA->At(di));
      os << "\tpdgId: " << dgp->pdgId << ", status: " << dgp->status << endl; 
    }
#endif
  }
}
void AnaBase::readPileUpHist() {
  size_t found = _puHistFile.find(".root");
  if (found == string::npos) {
    cerr << ">>> Warning: <<" << _puHistFile << ">> does not have .root extension!!" << endl;
    return;
  }

  const char* fname = gSystem->ExpandPathName(_puHistFile.c_str());
  if (gSystem->AccessPathName(fname)) {
    cerr << ">>> Warning: File <<" << _puHistFile << ">> was not found!!" << endl;
    return;
  }

  TFile file(fname);
  TH1D *h = dynamic_cast<TH1D*>(file.Get("plot_data_div_MC"));
  int nx = h->GetXaxis()->GetNbins();
  for (int i = 0; i < nx; ++i) {
    double wt = h->GetBinContent(i);
    _puWtList.push_back(wt);
  }
}
double AnaBase::wtPileUp(int& nPU) const {
  nPU = 0;

  const Event* evt = dynamic_cast<Event*>(eventA->At(0));
  assert(evt);
  std::vector<int> list = evt->nPU;
  if (!list.size()) return 1.0;

  int nbins = _puWtList.size();
  nPU = list.at(0);
  if (nPU < 0) nPU = 0;
  if (nPU >= nbins) nPU = nbins - 1;
  return _puWtList.at(nPU);
}
bool AnaBase::isTriggered(bool verbose) const {
  bool flag = false;
  for (size_t i = 0; i < _hltpaths->size(); ++i){
    string str = (*_hltpaths).at(i);
    bool found = false;
    //vector<string>::const_iterator it = find(_trigPathList.begin(), _trigPathList.end(), str);
    for (vector<string>::const_iterator it  = _trigPathList.begin();
	                                it != _trigPathList.end(); ++it) {
      if (str.find(*it) != string::npos) {
	found = true;
        break;
      }
    }
    if (!found) continue;

    if (verbose) cout << ">>> HLT Path = " << str 
                      << ", fired=" << (*_hltresults).at(i) 
                      << ", prescale=" << (*_hltprescales).at(i) << endl;
    if ((*_hltresults).at(i) == 1 && (*_hltprescales).at(i) == 1) {
      flag = true;
      break;
    }
  }
  return flag;
}
