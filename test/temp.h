#include "../atom.h"
#include "../vecr.h"
#include "../molecule.h"


template <class N>
class Number {

	private:
		N _number;

	public:
		Number (N i) {
			_number = i;
		}
}






Number<double> n (4.5);

