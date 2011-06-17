
import FWCore.ParameterSet.Config as cms
import copy

from PhysicsTools.PatAlgos.tools.coreTools import *
import PhysicsTools.PatAlgos.tools.helpers as patutils


def _setattr_ifexists(obj, attrName, attrValue):
	if hasattr(obj, attrName):
		setattr(obj, attrName, attrValue)

def switchToData(process):

#remove MC matching from standard PAT sequences
#from VHTauTau.TreeMaker.SwitchToData import switchToData
# remove MC matching

  removeMCMatching(process, ["All"], outputInProcess = False)
  removeMCMatching(process, ['METs'], "TC", outputInProcess = False)
  removeMCMatching(process, ['METs'], "PF", outputInProcess = False)

  process.patDefaultSequence.remove(process.patJetPartonMatch)
  process.patDefaultSequence.remove(process.patJetPartonMatchAK5PF)
  process.patDefaultSequence.remove(process.patJetGenJetMatchAK5PF)	
  process.patDefaultSequence.remove(process.patJetFlavourId)
  process.patDefaultSequence.remove(process.patJetPartons)
  process.patDefaultSequence.remove(process.patJetPartonAssociation)
  process.patDefaultSequence.remove(process.patJetPartonAssociationAK5PF)
  process.patDefaultSequence.remove(process.patJetFlavourAssociation)
  process.patDefaultSequence.remove(process.patJetFlavourAssociationAK5PF)
  
#runOnData(process, ["Jets"], outputInProcess = False)      

