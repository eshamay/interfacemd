#ifndef VECR_H_
#define VECR_H_

#include <iostream>
#include <math.h>
#include "matrixr.h"
#include "utility.h"

enum coord {x=0, y=1, z=2};

class MatR;

class VecR {

friend class MatR;

protected:
	std::vector<double>		_coords;

public:
	VecR ();
	VecR (const double x, const double y, const double z);
	VecR (const VecR& oldVec);						// A copy constructor
	VecR (const std::vector<double>& oldVec);
	VecR (const double * vec);
	~VecR ();

// Operators
	VecR 	operator+ (const VecR& input) const;						// Vector addition
	VecR 	operator- (const VecR& input) const;					// Vector subtraction
	void 	operator+= (const VecR& input);			// vector addition (assignment)
	void 	operator-= (const VecR& input);			// vector subtraction (assignment)
	double 	operator* (const VecR& input) const;		// Vector inner-product (dot-product)
	VecR  	operator* (const double input) const;	// Vector scaling
	VecR 	operator* (const MatR& input) const;
	void 	operator*= (const double input);		// Vector scaling (assignment)
	VecR	operator% (const VecR& input) const;		// Vector cross-product
	double	operator< (const VecR& input) const;		// Find the cos(angle) between two vectors
	double	operator[] (const coord index) const;	// Return the coordinate
	double	operator[] (const int index) const;		// Return the coordinate
	//void	operator= (VecR input);					// Set one vector equal to another
	bool	operator== (const VecR& input) const;	// Check identity

// Input & vector manipulation
	void Set (double X, double Y, double Z) {
		//printf ("_coords in (%d)\n", _coords);
		_coords[0] = X; // set all elements of a vector
		_coords[1] = Y;
		_coords[2] = Z;
	}
	void Set (const coord axis, const double val) { _coords[axis] = val; }

	void Zero ();								// Zero all elements of a vector
	void Scale (double val);					// scale the entire vector's magnitude
	void Scale (VecR val);						// scale each individual element by each element of another vector
	VecR MinVector (const VecR& input, const VecR& size) const;	// Find the min offset vector between two points in periodic bounds
	double MinDistance (const VecR& input, const VecR& size) const;
	//VecR RotateToFrame (VecR const * const frame) const; 
	VecR Wrap (VecR size, VecR origin = VecR ());						// Used to wrap a vector into a central periodic cell of the given size

	void X (const double val) { _coords[x] = val; }
	void Y (const double val) { _coords[y] = val; }
	void Z (const double val) { _coords[z] = val; }

// Output
	double Magnitude () const;
	VecR Unit () const;						// returns a unit vector in the same direction of this vector
	void Print () const;
};
#endif
