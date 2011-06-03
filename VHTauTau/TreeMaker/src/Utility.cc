#include <cassert>
#include "TTree.h"
#include "TROOT.h"
#include "VHTauTau/TreeMaker/interface/Utility.h"

TTree* Utility::getTree(const std::string& tree_name) {
  TTree* tree = dynamic_cast<TTree*>(gROOT->FindObject(tree_name.c_str()));
  assert(tree);
  return tree;  
}
