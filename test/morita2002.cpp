#include "morita2002.h"

using namespace std;

// used to calculate the SFG spectrum based on the morita/hynes 2002 method
int main (int argc, char **argv) {

	WaterSystemParams params ("morita2002.dat", 25000);

	//Morita2002<XYZSystem> sfg (params, "pos.xyz", VecR(12.0,12.0,20.0), "wanniers");
	Morita2002<AmberSystem> sfg (params, "prmtop", "mdcrd", "mdvel");

	sfg.Analyze();

return 0;
}
