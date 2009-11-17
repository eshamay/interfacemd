#ifndef UTILITY_H_
#define UTILITY_H_
#include <vector>
#include <algorithm>

// These RUN macros should clean up code for checking through vectors and arrays
#ifndef RUN
#define RUN(list)   \
  for (size_t i = 0; i < list.size(); i++)
#endif

#ifndef RUN2
#define RUN2(list)  \
  for (size_t j = 0; j < list.size(); j++)
#endif

#ifndef REMOVE_IF
#define REMOVE_IF(vec,pred)	\
  vec.erase (std::remove_if (vec.begin(), vec.end(), pred()), vec.end() );
#endif

// A general method for making new functors
#ifndef MAKE_FUNCTOR
#define MAKE_FUNCTOR(name,return_type,arg_type,code)	\
  class name {	\
	public:	\
	  return_type operator() (arg_type t) { code }	\
  };
#endif

// Creates a new predicate object that has a name as given, and takes a value of the given type (can be an object), and returns the value of the 'code' provided. code must return a boolean for this to work.
#ifndef MAKE_PREDICATE
#define MAKE_PREDICATE(name,arg_type,code)	\
  MAKE_FUNCTOR(name, bool, arg_type, return code;);
#endif

#ifndef FOR_EACH
#define FOR_EACH(name,fn)	\
  std::for_each (name.begin(), name.end(), fn);
#endif

// Maps OP onto A and stores the result into B
#ifndef MAP_TO
#define MAP_TO(A,B,OP)	\
  std::transform(A.begin(), A.end(), B.begin(), OP);
#endif

// Maps OP onto A and stores the result back into A - destructive function
#ifndef MAP
#define MAP(A,OP)	\
  MAP_TO(A,A,OP);
#endif
#endif
