The code and configuration needed to create the ntuple for Phase2 Level1 EGamma studies are kept here. The ntuple basically stores Phase2 Lavel1 objects namely

## Directory Structure of TTAnalysis
It can easily pulled within the L1Trigger package of CMSSW and contains

* interface/src : definition of objects stored in the ntuple
* python : configuration of the ntuple maker
* test : the ntuple maker TTAnalysis.cc 

The collection stored are

* Emulated L1TTTracks
* L1Egamma candidates
* L1TkElectron and L1EM (photon) candidates
* L1 PF and L1Puppi candidates
* L1 EmulatedVertex candidates

one can pull teh package inside L1Trigger subsystem of CMSSw to run it.

