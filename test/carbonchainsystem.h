#ifndef	CARBONCHAINSYSTEM_H_
#define	CARBONCHAINSYSTEM_H_

#include "../analysis.h"
#include "../carbonchain.h"
#include "../decane.h"
#include "../utility.h"

class CarbonChainSystem : public Analyzer {

  public:

	CarbonChainSystem 
	  (
	   const void * analysis_params,
	   const WaterSystemParams& params
	  );
	virtual ~CarbonChainSystem ();

};

#endif
