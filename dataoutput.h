#ifndef DATAOUTPUT_H_
#define DATAOUTPUT_H_

#include "patterns.h"
#include <iostream>

namespace md_analysis {

	class StarStatusBarUpdater : public patterns::observer::observer {
		protected:
			virtual void _updateScreenStatus ();
			long int _frequency;
			long int _count;
			long int _maxcount;

		public:
			StarStatusBarUpdater () : _frequency(0), _count(0), _maxcount(0) { }
			virtual ~StarStatusBarUpdater () { }

			StarStatusBarUpdater (const int frequency, const int maxcount, const int startingcount = 0) 
				: _frequency(frequency), _count(startingcount), _maxcount(maxcount) { }

			void Set (const int frequency, const int maxcount, const int startingcount = 0) {
				_frequency = frequency;
				_maxcount = maxcount;
				_count = startingcount;
			}

			// every time the updater is called the count is updated, and then output is performed based on the specific frequency supplied
			virtual void notify () {
				_count++;
				this->_updateScreenStatus ();
			}

	};	// Star status bar updater

}	// namespace md_analysis


#endif
