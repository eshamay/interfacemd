#ifndef PATTERNS_H_
#define PATTERNS_H_

#include <list>

namespace patterns {


	namespace observer {

		class observer {
			public:
				virtual void notify () = 0;
		};	// observer


		class observable {
			protected:
				typedef std::list<observer *> observer_list;
				observer_list observers;

			public:
				virtual void registerObserver (observer * o) {
					observers.push_back(o);
				}

				virtual void unregisterObserver (observer * o) {
					observers.remove(o);
				}

				virtual void notifyObservers () {
					for (observer_list::iterator it = observers.begin(); it != observers.end(); it++) {
						(*it)->notify();
					}
				}
		};	// observable


	}	// namespace observer



}	// patterns

#endif
