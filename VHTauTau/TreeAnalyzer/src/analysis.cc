#include <iostream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TStopwatch.h"

#include "AnaBase.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char* argv[]) {
   // Create  analysis object 
   AnaBase myana;

   // Read job input
   vector<std::string> fileList;
   fileList.push_back("DYToMuMu_1_2_T12.root");
   for (unsigned int i = 0; i < fileList.size(); i++)
     myana.setInputFile(fileList[i].c_str());

   gROOT->SetBatch(kTRUE);

   // Now go
   TStopwatch timer;
   cout << "Start event loop now with " << fileList.size() << " files" << endl;
   timer.Start();
   myana.eventLoop();
   timer.Stop();
   cout << "Realtime/CpuTime = " << timer.RealTime() << "/" << timer.CpuTime() << endl;

   return 0;
}
