import FWCore.ParameterSet.Config as cms

from VHTauTau.TreeMaker.EventBlock_cfi import eventBlock
from VHTauTau.TreeMaker.VertexBlock_cfi import vertexBlock
from VHTauTau.TreeMaker.JetBlock_cfi import jetBlock
from VHTauTau.TreeMaker.ElectronBlock_cfi import electronBlock
from VHTauTau.TreeMaker.METBlock_cfi import metBlock
from VHTauTau.TreeMaker.MuonBlock_cfi import muonBlock
from VHTauTau.TreeMaker.TauBlock_cfi import tauBlock
from VHTauTau.TreeMaker.GenParticleBlock_cfi import genParticleBlock
from VHTauTau.TreeMaker.GenJetBlock_cfi import genJetBlock
from VHTauTau.TreeMaker.GenMETBlock_cfi import genMETBlock
from VHTauTau.TreeMaker.TrackBlock_cfi import trackBlock
from VHTauTau.TreeMaker.TriggerBlock_cfi import triggerBlock
from VHTauTau.TreeMaker.PhotonBlock_cfi import photonBlock
from VHTauTau.TreeMaker.CommonVertexBlock_cfi import commonVertexBlock

treeContentSequence = cms.Sequence(
   eventBlock
 + vertexBlock
 + commonVertexBlock
 + electronBlock
 + genParticleBlock
 + genJetBlock
 + genMETBlock  
 + jetBlock  
 + metBlock
 + muonBlock
 + tauBlock
 + trackBlock
 + triggerBlock
 + photonBlock
)
