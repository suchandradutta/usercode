#ifndef __ANABASE__HH
#define __ANABASE__HH

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH1F.h"

#include "PhysicsObjects.h"

// REQUIRED, most probably to solve circular dependence problem!!!
class AnaBase;
class TClonesArray;
class TChain;
class TFile;


class VertexComparator {
public:
  bool operator()(const vhtm::Vertex &a, const vhtm::Vertex &b) const {
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

typedef struct  
{
  bool verbose;
  bool usesbit;
} Options;

class AnaBase {
    
public:

  AnaBase();
  virtual ~AnaBase();
    
  virtual void eventLoop() = 0;  // the main analysis 
  virtual bool beginJob();
  virtual void endJob() = 0;

  virtual void selectEvent() = 0;

  int setInputFile(const std::string& fname);
  int getEntries() const;
  void setData(int val);  
  int getRunNumber() const;
  virtual bool openFiles();
  virtual void closeFiles(); 
  virtual bool readJob(const std::string& jobFile, int& nFiles);
  virtual void printJob(std::ostream& os=std::cout) const;

  bool isTriggered(bool verbose=false) const;
  double wtPileUp(int& nPU) const;
  void readPileUpHist();
  
  void clearEvent();
  void enableBranches();
  int getEntry(int lflag) const;
  void setAddresses(); 

  void dumpGenInfo(std::ostream& os=std::cout) const;

  void findVtxInfo(std::vector<vhtm::Vertex>& list, std::map<std::string, double> cutMap, Options& op, std::ostream& os=std::cout) const;
  void findElectronInfo(std::vector<vhtm::Electron>& list, std::map<std::string, double> cutMap, Options& op, std::ostream& os=std::cout) const;
  void findMuonInfo(std::vector<vhtm::Muon>& list, std::map<std::string, double> cutMap, Options& op, std::ostream& os=std::cout) const;
  void findTauInfo(std::vector<vhtm::Tau>& list, std::map<std::string, double> cutMap, double vz, Options& op, std::ostream& os=std::cout) const;
  void findJetInfo(std::vector<vhtm::Jet>& list, std::map<std::string, double> cutMap, Options& op, std::ostream& os=std::cout) const;
  void dumpEvent(const std::vector<std::map<std::string, double> >& maps, std::ostream& os=std::cout) const;

protected:
  TChain* _chain;                // chain contains a list of root files containing the same tree
   TFile* _histf;                // The output file with histograms

public:
  // The tree branches
  TClonesArray* eventA;
  TClonesArray* vertexA;
  TClonesArray* genEventA;
  TClonesArray* tauA;
  TClonesArray* electronA;
  TClonesArray* muonA;
  TClonesArray* jetA;
  TClonesArray* metA;
  TClonesArray* genParticleA;
  TClonesArray* genJetA;
  TClonesArray* genMetA;
  TClonesArray* triggerobjA;
  TClonesArray* trackA;

 // Number of objects in each event (and each TClonesArray)
  int n_vertex;
  int n_tau;
  int n_electron;
  int n_muon;
  int n_jet;
  int n_met;  
  int n_genjet;
  int n_genparticle;
  int n_genmet;
  int n_triggerobj;
  int n_track;

  std::vector<int>* _l1physbits;
  std::vector<int>* _l1techbits;
  std::vector<std::string>* _hltpaths;
  std::vector<int>* _hltresults;
  std::vector<int>* _hltprescales;

  std::vector<std::string> _brList;
  std::vector<double> _puWtList;
  double _puevWt;

public:
  int nEvents;

protected:
  ofstream _fLog;   
  ofstream _evLog;   

public:
  bool _isMC;
  bool _readTrk;
  bool _readTrigObject;
  std::vector<std::string> _fileList;
  int _logOption;
  bool _studyTrigger;
  std::vector<std::string> _trigPathList;

  std::string _histFile;
  std::string _puHistFile;
  std::string _logFile;
  std::string _evFile;
  int _maxEvt;
};
#endif
