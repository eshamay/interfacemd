#include <iostream>
#include "vecr.h"

/*
	int
main(void)
{
	double          *a = new double[4];
	a[0] = 0., a[1] = 1., a[2] = 2., a[3] = 2.;

	double          *b = new double[4];
	b[0] = 0.5, b[1] = 0.5, b[2] = 0.5, b[3] = 0.5;

	double          *c = new double[4];

	Eigen::Map < VecR > am(a);
	Eigen::Map < VecR > bm(b);
	Eigen::Map < VecR > cm(c);
	cm = am + bm;

	std::cout << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << std::endl;
	cm.Print();

	c[2] = 7.0;
	c[0] = 3.0;

	std::cout << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << std::endl;
	cm.Print();

	delete [] a;
	delete [] b;
	delete [] c;

	return 0;
}
*/




#include "test.h"
int main () {

	double * a = new double[9];
	for (int i = 0; i < 9; i++)
		a[i] = i*.5;

	/*********** save this here ***************/
	std::vector< Eigen::Map<VecR> > vm;
	for (int i = 0; i < 3; i++)
		vm.push_back(Eigen::Map<VecR> (&a[3*i]));
	
	for (int i = 0; i < 9; i++)
		printf ("%.1f ", a[i]);
	printf ("\n");
	for (int i = 0; i < 3; i++)
		vm[i].Print();

	a[1] = 7.0;
	a[5] = 7.0;
	a[8] = 8.0;


	for (int i = 0; i < 9; i++)
		printf ("%.1f ", a[i]);
	printf ("\n");
	for (int i = 0; i < 3; i++)
		vm[i].Print();

	delete[] a;

	return 0;

}
