#include "../utility.h"
#include "../analysis.h"

class GMXAnalyzer : public Analyzer<GMXSystem> {

  private: 
  public:
    GMXAnalyzer (WaterSystemParams& wsp) : Analyzer<GMXSystem>(wsp) { return; }

};
