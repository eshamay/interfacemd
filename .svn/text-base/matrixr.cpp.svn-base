#include "matrixr.h"

MatR MatR::operator+ (const MatR& input) const {
	return MatR (_elements + input._elements);
}

VecR MatR::operator* (const VecR& input) const {		// Vector rotation/matrix-vector inner product

	boost::numeric::ublas::vector<double> vec (3);

	for (int i = 0; i < 3; i++) {
		vec(i) = input[i];
	}

return VecR (prod(_elements, vec));
}
/*
	double a,b,c;

	VecR output;

	for (int i=0; i<3; i++) {
		double val = 0.0;
		for (int j=0; j<3; j++) {
			val += this->Index(i,j) * input[j];
		}
		output.Set (coord(i), val);
	}
	//a = _elements[0]*input[x] + _elements[3]*input[y] + _elements[6]*input[z];
	//b = _elements[1]*input[x] + _elements[4]*input[y] + _elements[7]*input[z];
	//c = _elements[2]*input[x] + _elements[5]*input[y] + _elements[8]*input[z];
}
*/

MatR MatR::operator* (const MatR& input) const {		// Matrix rotation/multiplication
/*
	MatR out;

	double val = 0.0;
	//double const * elements = input.Elements();

	for (int i = 0; i < 3; i++) {
		for (int k = 0; k < 3; k++) {
	
			val = 0.0;

			for (int j=0; j<3; j++)
				val += this->Index(i,j) * input.Index(j,k);

			out.Set (i, k, val);
		}
	}
*/
return (MatR (prod(_elements, input._elements)));
}

void MatR::Print () const {
	for (int row=0; row < 3; row++) {
		for (int col=0; col<3; col++) 
			printf ("% 8.4f\t", _elements(row, col));
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
/*
	MatR out;

	for (int row=0; row<3; row++) {
		for (int col=0; col<3; col++)  
			out.Set(row, col, _elements[row*3+col]);
	}
*/
return (MatR (trans(_elements)));
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
	MatR out = U.Inverse() * (*this) * U;

	return (out);
}
#endif

double MatR::Trace () const {
	
	double out = _elements(0,0) + _elements(1,1) + _elements(2,2);

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
