#ifndef VECR_H_
#define VECR_H_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <iostream>
#include <math.h>
//#include "matrixr.h"
#include "utility.h"
//#include "FTensor.h"

enum coord {x=0, y=1, z=2};

typedef std::vector<double>	Double_vector;
typedef Double_vector::iterator Double_it;

typedef boost::numeric::ublas::vector<int> int_vector_t;
typedef boost::numeric::ublas::vector<double> vector_t;

class VecR : public vector_t {

public:

	typedef double T;
	typedef boost::numeric::ublas::vector_range<vector_t> range_t;

	VecR ();
	VecR (const T X, const T Y, const T Z);
	VecR (const VecR& oldVec);						// A copy constructor
	VecR (const T * vec);
	~VecR ();


// Operators
	//VecR 	operator+ (const VecR& input) const;						// Vector addition
	VecR 	operator+ (const T input) const;
	VecR 	operator+ (const VecR& input) const;
	VecR 	operator- (const VecR& input) const;					// Vector subtraction
	VecR 	operator- (const double input) const;
	void 	operator+= (const VecR& input);			// vector addition (assignment)
	void 	operator+= (const double input);			// vector addition (assignment)
	void 	operator-= (const VecR& input);			// vector subtraction (assignment)
	void 	operator-= (const double input);
	T 		operator* (const VecR& input) const;		// Vector inner-product (dot-product)
	VecR  	operator* (const double input) const;	// Vector scaling
	//void 	operator*= (const double input);		// Vector scaling (assignment)
	//VecR 	operator/ (const double input) const;
	//void 	operator/= (const double input);
	VecR	operator% (const VecR& input) const;		// Vector cross-product
	double	operator< (const VecR& input) const;		// Find the cos(angle) between two vectors
	T	operator[] (const coord index) const { return (*this)(index); }
	T	operator[] (const int index) const { return (*this)(index); }
	//double operator() (const coord index) const;
	//double operator() (const int index) const;
	//VecR&	operator= (const VecR& input);					// Set one vector equal to another
	//bool	operator== (const VecR& input) const;	// Check identity

// Input & vector manipulation
	void Set (T X, T Y, T Z) {
	  (*this)(0) = X;
	  (*this)(1) = Y;
	  (*this)(2) = Z;
	}
	void Set (const coord axis, const T val) { (*this)(axis) = val; }

	void Zero ();								// Zero all elements of a vector

	void X (const T val) { (*this)(0) = val; }
	void Y (const T val) { (*this)(1) = val; }
	void Z (const T val) { (*this)(2) = val; }

// Output
	T Magnitude () const;
	VecR Unit () const;						// returns a unit vector in the same direction of this vector
	void Print () const;
};

//typedef std::vector<VecR> VecR_vec;
//typedef std::vector<VecR>::const_iterator VecR_it;

#endif
