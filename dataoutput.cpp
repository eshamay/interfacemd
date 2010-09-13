#include "dataoutput.h"

	void md_analysis::StarStatusBarUpdater::_updateScreenStatus () {
		if (!(_count % (this->_frequency * 10)))
			std::cout << std::endl << _count << "/" << this->_maxcount << " ) ";
		if (!(_count % this->_frequency)) {
			std::cout << "*";
		}

		fflush (stdout);
	}
