#include "matrixr.h"

MatR MatR::operator+ (const MatR& input) const {
	MatR m (_elements);
	for (unsigned int i = 0; i < 9; i++)
		m._elements[i] += input[i];
	return (m);
}

VecR MatR::operator* (const VecR& input) const {		// Vector rotation/matrix-vector inner product
	VecR v;

	for (unsigned int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < 3; j++) {
			v._coords[i] += (Index(i,j) * input[j]);
		}
	}
	return (v);
}

MatR MatR::operator* (const MatR& input) const {		// Matrix rotation/multiplication
	MatR m;
	m.Zero();

	double val;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {

			val = 0.0;
			for (int k = 0; k < 3; k++) {
				val += (Index(i,k) * input.Index(k,j));
			}
			m.Set(i,j,val);
		}
	}

return (m);
}

double	MatR::operator[] (element const index) const {	// Return the coordinate
	if (index < 0 || index > 8) {
		std::cout << "Trying to access illegal index in MatrixR::operator[] (int const index)\nFix This\nMatrix:" << std::endl;
		this->Print();
		exit(1);
	}
	return _elements[index];
}	

double	MatR::operator[] (int const index) const {	// Return the coordinate
	if (index > 8) {
		std::cout << "Trying to access illegal index in MatrixR::operator[] (int const index)\nFix This\nMatrix:" << std::endl;
		this->Print();
		exit(1);
	}
	return _elements[index];
}

double	MatR::Index (int const row, int const col) const {	// Return the element
	if (row > 2 || col > 2) {
		std::cout << "Trying to access illegal index in MatrixR::Index (int const row, int const col)\nFix This\nMatrix:" << std::endl;
		this->Print();
		exit(1);
	}
	return (_elements[row+3*col]);
}

void	MatR::Set (int const row, int const col, double const val) {	// Set the element

	if (row > 2 || row < 0 || col > 2 || col < 0) {
		std::cout << "Trying to access illegal index in MatrixR::Set (int const row, int const col, double const val)\nFix This\nMatrix:" << std::endl;
		this->Print();
		exit(1);
	}

	_elements[row+3*col] = val;

	return;
}

void MatR::Print () const {
	for (int row=0; row < 3; row++) {
		for (int col=0; col<3; col++) 
			printf ("% 8.4f\t", this->Index(row, col));
		printf("\n");
	}
}

#ifdef _LINALG_
void MatR::CalcEigenSystem () {
	
	double A[9];
	for (int i=0; i<9; i++)
		A[i] = _elements[i];

	char jobvl = 'N';
	char jobvr = 'V';
	int N = 3;
	double vl[9];
	int ldvl = 3;
	int ldvr = 3;
	double work[100];
	int lwork = 50;
	int lda = 3;
	int info;

	dgeev_(&jobvl, &jobvr, &N, A, &lda, _eigenvalsR, _eigenvalsI, vl, &ldvl, _eigenvecs, &ldvr, work, &lwork, &info);
	
	_eigenset = true;
}

vector<VecR> MatR::EigenVectors () {

	if (!_eigenset) this->CalcEigenSystem();

	vector<VecR> out;
	out.push_back (VecR(_eigenvecs[0], _eigenvecs[1], _eigenvecs[2]));
	out.push_back (VecR(_eigenvecs[3], _eigenvecs[4], _eigenvecs[5]));
	out.push_back (VecR(_eigenvecs[6], _eigenvecs[7], _eigenvecs[8]));
	
	return (out);
}

vector< complex<double> > MatR::EigenValues () {
	
	if (!_eigenset) this->CalcEigenSystem();

	vector< complex<double> > out;
	out.push_back (complex<double> (_eigenvalsR[0], _eigenvalsI[0]));
	out.push_back (complex<double> (_eigenvalsR[1], _eigenvalsI[1]));
	out.push_back (complex<double> (_eigenvalsR[2], _eigenvalsI[2]));

	return (out);
}
#endif

MatR MatR::Transpose () const {
	MatR m;

	for (int row=0; row<3; row++) {
		for (int col=0; col<3; col++)  
			m.Set(row, col, this->Index(col,row));
	}
	return (m);
}

#ifdef _LINALG_
MatR MatR::Quaternion () {
	
	if (!_eigenset) this->CalcEigenSystem();

	MatR out (_eigenvecs);

	return(out);
}

MatR MatR::Inverse () const {
	
	int m = 3;
	int n = 3;
	int lda = 3;
	double A[9];
	double work[100];
	int lwork = 100;
	int info;
	int ipiv[9];

	// load up the elements of the array to be factored
	for (int i=0; i<9; i++)
		A[i] = _elements[i];

	// now call the L & U factorization routine from lapack
	dgetrf_ (&m, &n, A, &lda, ipiv, &info);

	dgetri_ (&n, A, &lda, ipiv, work, &lwork, &info);

	MatR out (A);

	return(out);
}

MatR MatR::Diagonalize () {
	
	if (!_eigenset) this->CalcEigenSystem();

	MatR U = this->Quaternion();
	MatR out = (U.Inverse() * (*this) * U);

	return (out);
}
#endif

double MatR::Trace () const {
	
	double out = (_elements[xx] + _elements[yy] + _elements[zz]);

	return(out);
}

/*
// This matrix is now rotated to another frame, thus we supply the x, y, and z axes vectors that we want to rotate to.
MatR MatR::RotateToFrame (VecR const * const frame) const {
	
	VecR _x = frame[0];
	VecR _y = frame[0];
	VecR _z = frame[0];

	// here's the lab-frame coordinates that we rotate from
	VecR X (1.0, 0.0, 0.0);
	VecR Y (0.0, 1.0, 0.0);
	VecR Z (0.0, 0.0, 1.0);

	// now we build our direction cosine matrix (each element is the cosine of the angle between two axes)

	double rotation[9] = { _x<X, _x<Y, _x<Z, _y<X, _y<Y, _y<Z, _z<X, _z<Y, _z<Z };

	MatR rot (rotation);

return (rot * (*this));
}
*/
