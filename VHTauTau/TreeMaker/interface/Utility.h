#ifndef __TreeMaker_Utility_hh
#define __TreeMaker_Utility_hh

#include <string>

class TTree;

class Utility 
{
public:
  static TTree* getTree(const std::string& tree_name="vhtree");
};
#endif
