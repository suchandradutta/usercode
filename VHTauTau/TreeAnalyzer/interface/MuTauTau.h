#ifndef __MUTAUTAU__HH
#define __MUTAUTAU__HH

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"

#include "PhysicsObjects.h"
#include "AnaBase.h"

using namespace vhtm;

class MuTauTau : public AnaBase {
    
public:

  MuTauTau();
  virtual ~MuTauTau();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void selectEvent();

  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout) const;
  
  void clearLists();
  void findGenInfo();
  
  //  void GenDetLevelMatch (const std::vector<Muon>& muoList, const std::vector<GenParticle>& genMuonList);
  void computeDeltaR(const std::vector<Muon>& muoList, const std::vector<Tau>& tauList);

  virtual void bookHistograms();

public:
  double nEvtSel[15];
  int count_sig_event;

  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::Muon> muoList;
  std::vector<vhtm::Electron> eleList;
  std::vector<vhtm::Tau> tauList;
  std::vector<vhtm::Jet> bjetList;

  std::vector<vhtm::GenParticle> genMuonList;
  std::vector<vhtm::GenParticle> genTauList;
  std::vector<vhtm::GenParticle> GEN_MUON_LIST;
  std::vector<vhtm::GenParticle> genMuoaList;
  std::vector<vhtm::GenParticle> genMuobList;
  std::vector<vhtm::GenParticle> genHList;
  std::vector<vhtm::GenParticle> genWList;
  std::vector<vhtm::GenParticle> genTau_H_List;
  std::vector<vhtm::GenParticle> genMuon_Z_List;

public:
  std::map<std::string, double> _vtxCutMap;
  std::map<std::string, double> _muonCutMap;
  std::map<std::string, double> _electronCutMap;
  std::map<std::string, double> _tauCutMap;
  std::map<std::string, double> _bjetCutMap;
  std::map<std::string, double> _evselCutMap;
};
#endif
