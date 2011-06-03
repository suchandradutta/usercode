#ifndef __CINT__
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TStyle.h>
#include <TTree.h>
#include <TLegend.h>
#include "VHAnalyser.h"
#endif
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

//int main() {
void create_histograms() {  
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include");
  gSystem->CompileMacro("VHAnalyser.C", "f");
  
  TFile* file = new TFile("MC_RelValTTbar.root"); 
  TTree* tree = (TTree*) file->Get("treeCreator/vhtree");
  
  std::string histo_fname = "histo_file.root"; 

  VHAnalyser* myanal = new VHAnalyser(tree, histo_fname);

  myanal->bookHistograms();
  myanal->Loop();
  myanal->saveHistograms();
  file->Close();
  delete myanal;
}

