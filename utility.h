#ifndef UTILITY_H_
#define UTILITY_H_
#include <vector>

// These RUN macros should clean up code for checking through vectors and arrays
#ifndef RUN

#define RUN(list)   \		// works for vectors - not for simple arrays
	for (size_t i = 0; i < list.size(); i++)

#endif

#ifndef RUN2

#define RUN2(list)  \
	for (size_t j = 0; j < list.size(); j++)

#endif

#define PRINTV(vector)	\
	printf ("vector size = % 5d\n", vector.size());

#define PRINTV2(vector)	\
	printf ("vector size = % 5d x % 5d\n", vector.size(), vector[0].size());

#endif
