#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>

int main () {
    using namespace boost::numeric::ublas;
    //slice s (0, 3, 3);
    //for (unsigned i = 0; i < s.size (); ++ i) {
        //std::cout << s (i) << " ";
    //}
	//std::cout << std::endl;
// definitions of matrices and vectors
    matrix<double> m (9,9);
    vector<double> v (3);

// setting/accessing matrix elements
	for (int i = 0; i < 9; i++) {
		//v[i] = (double)i;
	for (int j = 0; j < 9; j++) {
		m(i,j) = (i+1)*j;
	}
	}
    std::cout << v << std::endl;
    std::cout << m << std::endl;

// picking out sub-matrices from m
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			matrix_range<matrix<double> > mr (
				m,
				range (i*3, (i*3)+3),
				range (j*3, (j*3)+3)
			);
			std::cout << mr << std::endl;
		}
	}
return 0;
}
