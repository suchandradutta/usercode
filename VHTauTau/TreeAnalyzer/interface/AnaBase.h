#ifndef __ANABASE__HH
#define __ANABASE__HH

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include "TLorentzVector.h"
#include "TVector.h"

// REQUIRED, most probably to solve circular dependence problem!!!
class AnaBase;
class TClonesArray;
class TChain;
class TFile;

class AnaBase {
    
public:

  AnaBase(const std::string& filename="pippo.out");
  virtual ~AnaBase();
    
    void eventLoop();  // the main analysis 
    void setInputFile(const std::string& fname);
     int getEntries() const;
    void setData(int val);  
     int getRunNumber() const;
    bool openFiles();
    void closeFiles(); 
    bool readJob(const std::string& jobFile, int& nFiles);
    void printJob(std::vector<std::string>& fileList, std::ostream& os=std::cout);
  
  template <class T>
  bool fillHist1D(const std::string& hname,  T value, double w=1.0);
  template <class T1, class T2>
  bool fillHist2D(const std::string& hname,  T1 xvalue, T2 yvalue, double w=1.0);
  bool fillProfile(const std::string& hname, float xvalue, float yvalue, double w=1.0);

  virtual bool selectEvent();
  virtual void bookHistograms();
  void clearEvent();
  void enableBranches();
   int getEntry(int lflag) const;
  void setAddresses(); 

protected:
  void getBranches(int lflag) {};
  bool beginJob();
  void endJob();

private:
  TChain* _chain;                // chain contains a list of root files containing the same tree
   TFile* _histf;                // The output file with histograms

public:
  // The tree branches
  TClonesArray* eventA;
  TClonesArray* triggerA;
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

public:
  int nEvents;
  int nEvt[10];

protected:
  ofstream _fLog;   

public:
  int _logOption;
  std::map<std::string, double> _vtxCutMap;
  std::map<std::string, double> _muonCutMap;
  std::map<std::string, double> _electronCutMap;
  std::map<std::string, double> _tauCutMap;
  std::map<std::string, double> _bjetCutMap;
  std::string _histFile;
  std::string _logFile;
  int _maxEvt;
};
#endif
