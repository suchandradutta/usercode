#include "L1Trigger/TTAnalysis/interface/AnalObjects.h"
#include <iostream>

TTStudy::Event::Event() :
  event(-1),
  run(-1),
  nPileUp(-1),
  beamSpotX0(-999.9),
  beamSpotY0(-999.9),
  beamSpotZ0(-999.9),
  nL1Vtx(-1),
  zL1Vtx(-999.9),
  sumPtL1Vtx(-999.9),
  nL1EmuVtx(-1),
  zL1EmuVtx(-999.9),
  sumPtL1EmuVtx(-999.9),
  nOffVtx(-1),
  zOffVtx(-999.9),
  sumPtOffVtx(-999.9)
{}

// Information about SimTracks
TTStudy::SimTrack::SimTrack() :
  pt(-999.9),
  eta(-999.9),
  phi(-999.9),
  vx(-999.9),
  vy(-999.9),
  vz(-999.9),
  vtxIndx(-1),
  type(-1)
{}
// Information about Reconstruceted Tracks
TTStudy::Track::Track() :
  pt(-999.9),
  eta(-999.9),
  phi(-999.9),
  localPhi(-999.9),
  chiSquare(-999.9),
  chiSquareRed(-999.9),
  chi2RPhiRed(-999.9),
  chi2RZRed(-999.9),
  chi2Bend(-999.9),
  chi2BendRed(-999.9),
  curvature(-999.9),
  vertexX(-999.9),
  vertexY(-999.9),
  vertexZ(-999.9),
  nStub(-1),
  nStub_PS(-1),
  nStub_SS(-1),
  d0(-999.9),
  z0(-999.9),
  mva1(-999),
  mva2(-999),
  mva3(-999),
  nFitPars(-1),
  hitPattern(0)
{} 

// Information about Reconstruceted Tracks
TTStudy::OfflineTrack::OfflineTrack() :
  pt(-999.9),
  eta(-999.9),
  phi(-999.9),
  curvature(-999.9),
  chiSquare(-999.9),
  chiSquareRed(-999.9),
  vertexX(-999.9),
  vertexY(-999.9),
  vertexZ(-999.9),

  d0(-999.9),
  z0(-999.9),
  d0Err(-999.9),
  z0Err(-999.9),

  d0PV(-999.9),
  z0PV(-999.9),
  d0ErrPV(-999.9),
  z0ErrPV(-999.9)
{} 
TTStudy::GenParticle::GenParticle() :
  eta(-999.9),
  phi(-999.9),
  p(-999.9),
  px(-999.9),
  py(-999.9),
  pz(-999.9),
  pt(-999.9),
  energy(-999.9),
  vx(-999.9),
  vy(-999.9),
  vz(-999.9),
  status(-999),
  pdgId(-1),
  charge(-999),
  motherIndex(-999)
{
  daughterIndices.clear();
}  
TTStudy::L1Object::L1Object() :
  e(-999.9),
  et(-999.9),
  pt(-999.9),
  eta(-999.9),
  phi(-999.9),
  label(9999),
  hwQual(-999)
{}
TTStudy::L1TkObject::L1TkObject() :
  e(-999.9),
  et(-999.9),
  pt(-999.9),
  eta(-999.9),
  phi(-999.9),
  pt_tk(-999.9),
  eta_tk(-999.9),
  phi_tk(-999.9),
  curvature(-999.9),
  chiSquare(999.9),
  chiSquareRed(999.9),
  vertexX(-999.9),
  vertexY(-999.9),
  vertexZ(-999.9),

  d0(-999.9),
  z0(-999.9),    
  d0Err(-999.9),
  z0Err(-999.9),

  trackIsolation(-999.9),
  trackIsolationPV(-999.9),
  pfIsolation(-999.9),
  pfIsolationPV(-999.9),
  puppiIsolation(-999.9),
  puppiIsolationPV(-999.9),
  
  label(99),
  hwQual(-99)
{}
TTStudy::L1PFObject::L1PFObject() :
  e(-999.9),
  pt(-999.9),
  eta(-999.9),
  phi(-999.9),
  dxy(-999.9),
  z0(-999.9),
  puppyWeight(-999.0),
  pdgId(-999),    
  pt_refTrk(-999.0),
  eta_refTrk(-999.0),
  phi_refTrk(-999.0),
  pt_refClus(-999.0),
  eta_refClus(-999.0), 
  phi_refClus(-999.0)
{}
