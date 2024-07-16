#include "L1Trigger/TTAnalysis/interface/AnalObjects.h"
#include <vector>

namespace {
  struct dictionary {
    TTStudy::Event          rv1;
    TTStudy::SimTrack       rv2;
    TTStudy::Track          rv3;
    TTStudy::GenParticle    rv4;
    TTStudy::L1Object       rv5;
    TTStudy::L1TkObject     rv6;
    TTStudy::OfflineTrack   rv7;
    TTStudy::L1PFObject   rv8;
    std::vector<TTStudy::SimTrack>       vrv1;
    std::vector<TTStudy::Track>          vrv2;
    std::vector<TTStudy::GenParticle>    vrv3;
    std::vector<TTStudy::L1Object>       vrv4;
    std::vector<TTStudy::L1TkObject>     vrv5;
    std::vector<TTStudy::OfflineTrack>   vrv6;
    std::vector<TTStudy::L1PFObject>   vrv7;
  };
}
