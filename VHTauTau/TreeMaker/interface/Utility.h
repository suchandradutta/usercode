#ifndef Utility_hh
#define Utility_hh

#include <string>
class Tree;

class Utility 
{
public:
  static TTree* getTree(const std::string& tree_name="vhtree"); 
};
#endif
