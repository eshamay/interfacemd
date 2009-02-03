#ifndef UTILITY_H_
#define UTILITY_H_

using namespace std;

// These RUN macros should clean up code for checking through vectors and arrays
#ifndef RUN

#define RUN(list)   \
   	for (unsigned int i = 0; i < list.size(); i++)

#endif

#ifndef RUN2

#define RUN2(list)  \
   	for (unsigned int j = 0; j < list.size(); j++)

#endif

#define PRINTV(vector)	\
	printf ("vector size = % 5d\n", vector.size());

#define PRINTV2(vector)	\
	printf ("vector size = % 5d x % 5d\n", vector.size(), vector[0].size());

#endif

typedef std::vector<Atom>::iterator ATOM_IT;
typedef std::vector<Atom *>::iterator PATOM_IT;
typedef std::vector<Molecule *>::iterator PMOL_IT;
typedef std::vector<Molecule>::iterator MOL_IT;
typedef std::vector<double>::iterator DBL_IT;
