import FWCore.ParameterSet.Config as cms

from VHTauTau.TreeMaker.EventBlock_cfi import *
from VHTauTau.TreeMaker.VertexBlock_cfi import *
from VHTauTau.TreeMaker.CaloJetBlock_cfi import *
from VHTauTau.TreeMaker.JetBlock_cfi import *
from VHTauTau.TreeMaker.ElectronBlock_cfi import *
from VHTauTau.TreeMaker.METBlock_cfi import *
from VHTauTau.TreeMaker.MuonBlock_cfi import *
from VHTauTau.TreeMaker.TauBlock_cfi import *
from VHTauTau.TreeMaker.GenParticleBlock_cfi import *
from VHTauTau.TreeMaker.GenJetBlock_cfi import *
from VHTauTau.TreeMaker.GenMETBlock_cfi import *
from VHTauTau.TreeMaker.TriggerBlock_cfi import *
from VHTauTau.TreeMaker.SuperClusterBlock_cfi import *
from VHTauTau.TreeMaker.PhotonBlock_cfi import *

treeContentSequence = cms.Sequence(
   eventBlock
 + vertexBlock
 + caloJetBlock
 + jetBlock  
 + electronBlock
 + metBlock
 + muonBlock
 + tauBlock
 + genParticleBlock
 + genJetBlock
 + genMETBlock  
 + triggerBlock
 + photonBlock
)
