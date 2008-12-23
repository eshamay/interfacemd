#include "../matrixr.h"
#include "../vecr.h"

using namespace std;

int main (int argc, char*argv) {

	double a [9] = { 0,3,7,8,2,9,5,4,1 };
	double b [9] = { 9,2,1,8,7,5,3,4,6 };
	MatR ma;
	ma.Set(a);
	MatR mb (b);
	VecR va (6.,3.,2.);

	ma.Print();


return 0;
}
