#include "matrixr.h"

MatR MatR::operator+ (const MatR& input) const {
	MatR m;
	double val;
	for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
		val = _matrix[i][j] + input._matrix[i][j];
		m.Set(i,j,val);
	}}

	return (m);
}

// multiply a vector by a matrix
VecR MatR::operator* (const VecR& input) const {		// Vector rotation/matrix-vector inner product
	VecR v;
	for (int i = 0; i < 3; i++) {
		double val = 0.0;
		for (int j = 0; j < 3; j++) {
			val += _matrix[i][j] * input[j];
		}
		v._coords[i] = val;
	}

	return (v);
}

// multiply a matrix by a matrix
MatR MatR::operator* (const MatR& input) const {		// Matrix rotation/multiplication
	MatR m;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			double val = 0.0;
			for (int k = 0; k < 3; k++) {
				val += _matrix[i][k] * input.Index(k,j);
			}
			m.Set(i,j,val);
		}
	}

	return (m);
}

double MatR::operator() (int const row, int const col) const {
	return (_matrix[row][col]);
}

double MatR::operator() (coord const row, coord const col) const {
	return (_matrix[row][col]);
}
double MatR::Index (int const row, int const col) const {	// Return the element
	return (_matrix[row][col]);
}
double MatR::Index (coord const row, coord const col) const {	// Return the element
	return (_matrix[row][col]);
}

void MatR::Set (int const row, int const col, double const val) {	// Set the element
	_matrix[row][col] = val;
}

void MatR::Set (coord const row, coord const col, double const val) {	// Set the element
	_matrix[row][col] = val;
}

// set the matrix using a pre-built array of data
void MatR::Set (double * const data) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			_matrix[i][j] = data[i*3+j];
		}
	}
}

void MatR::Print () const {
	for (int row=0; row < 3; row++) {
		for (int col=0; col<3; col++)
			printf ("% 8.4f\t", _matrix[row][col]);
		printf("\n");
	}
}

#ifdef _LINALG_
/*
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
*/
#endif

MatR MatR::Transpose () const {
	MatR m;
	for (unsigned int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < 3; j++) {
			m.Set(i,j,_matrix[j][i]);
		}
	}
	return (m);
}

MatR MatR::Inverse () const {

	double G[9];

	// load up the elements of the array to be factored
	for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
		G[i*3+j] = _matrix[i][j];
	}}

	int n = 3;
	int ipiv[n];
	int info;
	/*---------Call Cell LAPACK library---------*/
	dgetrf_(&n, &n, G, &n, ipiv, &info);
	if( info != 0 )
		std::cout << "hey! dgetrf didn't work :( - matrix factorization is dubious" << std::endl;

	/*---------Query workspace-------*/
	// find out how much space we need to do this thing
	double workspace;
	int tmp=-1;
	int lwork;
	dgetri_(&n, G, &n, ipiv, &workspace, &tmp, &info);

	lwork = (int)workspace;
	double work[lwork];
	//work = (double * ) malloc (sizeof(double)*lwork);
	/*---------Call Cell LAPACK library---------*/
	dgetri_(&n, G, &n, ipiv, work, &lwork, &info);
	if (info != 0)
		std::cout << "matrix inverse died after 2nd dgetri_ call" << std::endl;

	// G now holds the inverse

	MatR out (G);

	return(out);
}

#ifdef _LINALG_
/*
MatR MatR::Quaternion () {

	if (!_eigenset) this->CalcEigenSystem();

	MatR out (_eigenvecs);

	return(out);
}
*/

/*
MatR MatR::Diagonalize () {

	if (!_eigenset) this->CalcEigenSystem();

	MatR U = this->Quaternion();
	MatR out = (U.Inverse() * (*this) * U);

	return (out);
}
*/
#endif

double MatR::Trace () const {

	double out = _matrix[0][0] + _matrix[1][1] + _matrix[2][2];

	return(out);
}

double MatR::Determinant () const {

	double det = 0.0;
	det += _matrix[0][0]*_matrix[1][1]*_matrix[2][2];
	det -= _matrix[0][0]*_matrix[1][2]*_matrix[2][1];
	det += _matrix[0][1]*_matrix[1][2]*_matrix[2][0];
	det -= _matrix[0][1]*_matrix[1][0]*_matrix[2][2];
	det += _matrix[0][2]*_matrix[1][0]*_matrix[2][1];
	det -= _matrix[0][2]*_matrix[1][1]*_matrix[2][0];

	return (det);
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
