// This is the dipole field tensor (T) from the Morita+Hynes paper from 2002 (JPC B, 106, 3, 2002, p.675)

class DipoleFieldTensor {

private:
	double const * _r;	// distance vector (x,y,z) between 2 molecules
	double _tensor[9];	// The internal elements of the dipole field tensor
	double _distance;	// scalar distance between the two molecules

public:

	DipoleFieldTensor (double const * const r);		// the dipole field tensor only depends on the distance vector
	~DipoleFieldTensor ();

	// elements are indexed as [row][column]
	double const Element (int const row, int const col) const { return _tensor[col*3+row]; }
	double const * Tensor () { return _tensor; }
};

