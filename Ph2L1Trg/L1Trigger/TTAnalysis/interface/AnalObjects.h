/// ////////////////////////////////////////
/// Objects to be stored in the ntuble
/// ////////////////////////////////////////

#ifndef ANALOBJECTS_H
#define ANALOBJECTS_H

#include "TObject.h"
#include <vector>

namespace TTStudy {
  class Event : public TObject {

  public:
    Event(); 
    ~Event(){}
    
    int event;
    int run;
    int nPileUp;
    float beamSpotX0;
    float beamSpotY0;
    float beamSpotZ0;
    int   nL1Vtx;
    float zL1Vtx;
    float sumPtL1Vtx;
    int   nL1EmuVtx;
    float zL1EmuVtx;
    float sumPtL1EmuVtx;
    int   nOffVtx;
    float zOffVtx;
    float sumPtOffVtx;

    ClassDef(Event, 1)
  };

  //  class Tracklet;
  class Tracklet : public TObject {
    
  public:
    Tracklet();
    ~Tracklet(){}

    int layer1;
    float pt1; 
    float phi1;
    float eta1;
    float x1;
    float y1;
    float z1;
    float r1;
    int trackIndex1;
    int particleId1;
    float truePt1;
    float deltaPhi1;
    float zIntercept1;
    
    int layer2;
    float pt2; 
    float phi2;
    float eta2;
    float x2;
    float y2;
    float z2;
    float r2;
    int trackIndex2;
    int particleId2;
    float truePt2;  
    float deltaPhi2;
    float zIntercept2;
    
    float phiMiss;
    float zMiss;
    float twoPointPt;
    float twoPointZIntercept;
    
    ClassDef(Tracklet, 1)
  };
  
  class SimTrack : public TObject {
    
  public:
    SimTrack();
    ~SimTrack(){}
    
    float pt;
    float eta;
    float phi;
    float vx;
    float vy;
    float vz;
    int vtxIndx;
    int type;
    
    ClassDef(SimTrack, 1)
  };

  class Track : public TObject {
    
  public:
    Track();
    ~Track() {}
    
    float pt;
    float eta;
    float phi;
    float localPhi;
    float chiSquare;
    float chiSquareRed;
    float chi2RPhiRed;
    float chi2RZRed;
    float chi2Bend;
    float chi2BendRed;
    float curvature;
    float vertexX;
    float vertexY;
    float vertexZ;

    int   nStub;
    int   nStub_PS;
    int   nStub_SS;

    float d0;
    float z0;    

    float mva1;
    float mva2;
    float mva3;
    
    int nFitPars;
    unsigned int hitPattern;

    ClassDef(Track, 1)
  };

  class OfflineTrack : public TObject {
    
  public:
    OfflineTrack();
    ~OfflineTrack() {}
    
    float pt;
    float eta;
    float phi;
    float curvature;
    float chiSquare;
    float chiSquareRed;
    float vertexX;
    float vertexY;
    float vertexZ;

    float d0;
    float z0;    
    float d0Err;
    float z0Err;

    float d0PV;
    float z0PV;
    float d0ErrPV;
    float z0ErrPV;

    ClassDef(OfflineTrack, 1)
  };

  class GenParticle : public TObject {
    
  public:
    GenParticle();
    ~GenParticle() {}
    
    float eta;
    float phi;
    float p;
    float px;
    float py;
    float pz;
    float pt;
    float energy;
    float vx;
    float vy;
    float vz;
    int status;
    int pdgId;
    int charge;
    int motherIndex;
    std::vector<int> daughterIndices;

    ClassDef(GenParticle, 1)
  };
  
  class L1Object : public TObject {
    
  public:
    L1Object();
    ~L1Object() {}
   
    float e;
    float et;
    float pt;
    float eta;
    float phi;
    unsigned int label;    
    int hwQual;
    std::vector<int> wpFlags;  
    ClassDef(L1Object, 1)
  };
  class L1TkObject : public TObject {
    
  public:
    L1TkObject();
    ~L1TkObject() {}
   
    float e;
    float et;
    float pt;
    float eta;
    float phi;
    float pt_tk;
    float eta_tk;
    float phi_tk;
    float curvature;
    float chiSquare;
    float chiSquareRed;
    float vertexX;
    float vertexY;
    float vertexZ;

    float d0;
    float z0;    
    float d0Err;
    float z0Err;
    
    float trackIsolation;
    float trackIsolationPV;
    float pfIsolation;
    float pfIsolationPV;
    float puppiIsolation;
    float puppiIsolationPV;

    unsigned int label;    
    int hwQual;

    ClassDef(L1TkObject, 1)
  };
  class L1PFObject : public TObject {
    
  public:
    L1PFObject();
    ~L1PFObject() {}
   
    float e;
    float pt;
    float eta;
    float phi;
    int charge;
    float dxy;
    float z0;    
    float puppyWeight;
    int   pdgId;
    float pt_refTrk;
    float eta_refTrk;
    float phi_refTrk;
    float pt_refClus;
    float eta_refClus;
    float phi_refClus;

    ClassDef(L1PFObject, 1)
  };
}
#endif
