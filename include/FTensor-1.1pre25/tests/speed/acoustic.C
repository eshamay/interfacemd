#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
using namespace std;

#include "../../FTensor.h"

void initial(double P1[], double P2[], double P3[],
	     double c[], const int N)
{
  int cr = N/2-1;
  int cc = 7.0*N/8.0-1;
  float s2 = 64.0 * 9.0 / pow(N/2.0,2.0);

  for(int i=0;i<N;++i)
    for(int j=0;j<N;++j)
      {
	c[i+N*j]=0.2;
	P1[i+N*j] = 0.0;
	P2[i+N*j] = exp(-(pow(i-cr,2.0)+pow(j-cc,2.0)) * s2);
	P3[i+N*j] = 0.0;
      }

  const int blockLeft = 0;
  const int blockRight = 2 * N / 5.0;
  const int blockTop = N / 3.0;
  const int blockBottom = 2 * N / 3.0;

  for(int i=blockTop; i<blockBottom; ++i)
    for(int j=blockLeft; j<blockRight; ++j)
      c[i+N*j] = 0.5;

  int channelLeft = 4*N/5.0;
  int channelRight = N;
  int channel1Height = 3*N/8.0;
  int channel2Height = 5*N/8.0;
  for(int i=channelLeft; i<channelRight; ++i)
    for(int j=channel1Height; j<channel2Height; ++j)
      c[i+N*j]=0;
}

void FTen(double P1_[], double P2_[], double P3_[],
	    double c_[], const int N)
{
  FTensor::Tensor0<double*> P1(P1_), P2(P2_), P3(P3_), c(c_);

  FTensor::Tensor1<int,2> d_ijk(1,N);
  FTensor::Tensor1<double,2> d_xyz(1,1);
  FTensor::Index<'l',2> l;
  FTensor::Index<'m',2> m;
  for(int j=0;j<N;++j)
    for(int i=0;i<N;++i)
      {
	if(i!=0 && i!=N-1 && j!=0 && j!=N-1)
	  {
//  	    P3 = 2*P2 - P1;
	    P3 = 2*P2 - P1
	      + c*(dd(P2,l,m,d_ijk,d_xyz)(0,0)
		   + dd(P2,l,m,d_ijk,d_xyz)(1,1));
	  }
	++P1;
	++P2;
	++P3;
	++c;
      }
//    cout << P3_[(N/2-1)+N*((7*N)/8-1)] << endl;
}


void simple(double P1[], double P2[], double P3[],
	    double c[], const int N)
{
  for(int j=1;j<N-1;++j)
    for(int i=1;i<N-1;++i)
      {
//  	P3[i+N*j] = 2*P2[i+N*j] - P1[i+N*j];
	P3[i+N*j] = (2-4*c[i+N*j])*P2[i+N*j] - P1[i+N*j]
	  + c[i+N*j]*(P2[i-1+N*j] + P2[i+1+N*j]
		      + P2[i+N*(j-1)] + P2[i+N*(j+1)]);
      }
//    cout << P3[(N/2-1)+N*((7*N)/8-1)] << endl;
}


int main(int argc, char *argv[])
{
  const int N=atoi(argv[1]);

  vector<double> P1(N*N), P2(N*N), P3(N*N), c(N*N);

  initial(&c[0],&P1[0],&P2[0],&P3[0],N);


  const int niters=atoi(argv[2]);

  for(int i=0;i<niters;++i)
    {
//        simple(&P1[0],&P2[0],&P3[0],&c[0],N);
//        simple(&P2[0],&P3[0],&P1[0],&c[0],N);
//        simple(&P3[0],&P1[0],&P2[0],&c[0],N);
      FTen(&P1[0],&P2[0],&P3[0],&c[0],N);
      FTen(&P2[0],&P3[0],&P1[0],&c[0],N);
      FTen(&P3[0],&P1[0],&P2[0],&c[0],N);
    }
}


