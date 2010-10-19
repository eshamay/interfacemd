#pragma once
#ifndef EIGEN_MATRIX_ADDON_H_
#define EIGEN_MATRIX_ADDON_H_

// for optimization
#ifndef NDEBUG
#define NDEBUG
#endif

typedef Eigen::MatrixBase<double> vector_t;
typedef double data_t;

// Arithmetic

template<typename OtherDerived>
inline OtherDerived	operator% (const OtherDerived& input) const		// Vector cross-product
{ return this->cross(input); }

template<typename OtherDerived>
inline double operator< (const OtherDerived& input) const		// Find the cos(angle) between two vectors
{
	// determine the cos(angle) between two vectors by applying the dot product, and dividing by the magnitudes
	// cos(angle) = dotproduct/magnitudes
	// Return cos(angle)
	return this->dot(input) / this->norm() / input.norm();
}

// Input & vector manipulation
inline void Set (const int axis, const data_t val) { this->operator[](axis) = val; }
//inline void	Set (int const row, int const col, double const val) { this->operator()(row,col) = val; }
inline void Set (const double X, const double Y, const double Z) { *this << X,Y,Z; }

// Output
inline data_t Magnitude () const { return this->norm(); }

template<typename OtherDerived>
inline OtherDerived Unit () const // returns the unit vector (normalized)
{ return this->normalized(); }

inline void Print () const
{
	if (cols() > 1)
		std::cout << *this << std::endl;
	else
		printf ("[% 8.3f% 8.3f% 8.3f ]\n", this->x(), this->y(), this->z());
}


template<typename OtherDerived>
OtherDerived& Transpose() const { return this->transpose(); }

#endif

