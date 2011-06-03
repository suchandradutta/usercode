#define VHAnalyser_cxx
#include <iostream>
#include "VHAnalyser.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void VHAnalyser::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L VHAnalyser.C
//      Root > VHAnalyser t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      std::cout << " Event # " << jentry << std::endl;
      std::cout << "nVertex " << nVertex << endl;
      std::cout << "nCaloJet " << nCaloJetrtex << endl;
      std::cout << "nJet " << nJet << endl;
      std::cout << "nElectron " << nElectron << endl;
      std::cout << "nMuon " << nMuon << endl;
      std::cout << "nTau " << nTau << endl;
      std::cout <<std::endl;
      // if (Cut(ientry) < 0) continue;
      hNVertex->Fill(nVertex);
   }
}
void VHAnalyser::bookHistograms() {
  if (!bookedHistos) {
    // Open Output File
    std::string root_file = fileName;
    outputFile = new TFile(root_file.c_str(), "RECREATE");
    outputFile->cd();
    hNVertex = new TH1F("nVertex", "Number of Vertex", 31, -0.5, 30.5);

  }
}
void VHAnalyser::saveHistograms() {
  if (outputFile) {
    hNVertex->Write();
    outputFile->Close();
  }
}
