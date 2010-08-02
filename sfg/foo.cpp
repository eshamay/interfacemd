#include "foo.h"

int main () {

  // for now, f is 'h', the 3Nx3 block identity tensor
  MatR A;
  A << 1,2,3,8,7,1,2,2,1;
  VecR _x;
  _x.setZero(3);
  VecR b (12,32,10);

  // now solve for f in the equation g*f = h 
  int N = 3;
  int nrhs = 1;
  int ipiv[N];
  for (int i = 0; i < N; i++) ipiv[i] = 0;

  double work[N*nrhs];
  float swork[N*(N+nrhs)];
  int iter;
  int info = 0;

  dsgesv (&N, &nrhs, &A(0,0), &N, ipiv, &b(0,0), &N, &_x(0,0), &N, work, swork, &iter, &info);

  A.Print();

  b.Print();

  _x.Print();

  A.lu().solve(b,&_x);

  A.Print();
  b.Print();
  _x.Print();

}
