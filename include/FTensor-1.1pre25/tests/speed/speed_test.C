/* Tests the difference that loop fusion makes. */

const int N=60000;

void fused(const double metric[6][N], double inverse[6][N])
{
  for(int i=0;i<N;++i)
    {
      double det=metric[0][i]*metric[3][i]*metric[5][i]
	+ metric[1][i]*metric[4][i]*metric[2][i]
	+ metric[2][i]*metric[1][i]*metric[4][i]
	- metric[0][i]*metric[4][i]*metric[4][i]
	- metric[1][i]*metric[1][i]*metric[5][i]
	- metric[2][i]*metric[3][i]*metric[2][i];
      inverse[0][i]=
	(metric[3][i]*metric[5][i] - metric[4][i]*metric[4][i])/det;
      inverse[1][i]=
	(metric[2][i]*metric[4][i] - metric[1][i]*metric[5][i])/det;
      inverse[2][i]=
	(metric[1][i]*metric[4][i] - metric[2][i]*metric[3][i])/det;
      inverse[3][i]=
	(metric[0][i]*metric[5][i] - metric[2][i]*metric[2][i])/det;
      inverse[4][i]=
	(metric[2][i]*metric[1][i] - metric[0][i]*metric[4][i])/det;
      inverse[5][i]=
	(metric[3][i]*metric[0][i] - metric[1][i]*metric[1][i])/det;
    }
}

void fused_big(const double metric[6][N], double inverse[6][N])
{
  double det[N];
  for(int i=0;i<N;++i)
    {
      det[i]=metric[0][i]*metric[3][i]*metric[5][i]
	+ metric[1][i]*metric[4][i]*metric[2][i]
	+ metric[2][i]*metric[1][i]*metric[4][i]
	- metric[0][i]*metric[4][i]*metric[4][i]
	- metric[1][i]*metric[1][i]*metric[5][i]
	- metric[2][i]*metric[3][i]*metric[2][i];
    }
  for(int i=0;i<N;++i)
    {
      inverse[0][i]=
	(metric[3][i]*metric[5][i] - metric[4][i]*metric[4][i])/det[i];
    }
  for(int i=0;i<N;++i)
    {
      inverse[1][i]=
	(metric[2][i]*metric[4][i] - metric[1][i]*metric[5][i])/det[i];
    }
  for(int i=0;i<N;++i)
    {
      inverse[2][i]=
	(metric[1][i]*metric[4][i] - metric[2][i]*metric[3][i])/det[i];
    }
  for(int i=0;i<N;++i)
    {
      inverse[3][i]=
	(metric[0][i]*metric[5][i] - metric[2][i]*metric[2][i])/det[i];
    }
  for(int i=0;i<N;++i)
    {
      inverse[4][i]=
	(metric[2][i]*metric[1][i] - metric[0][i]*metric[4][i])/det[i];
    }
  for(int i=0;i<N;++i)
    {
      inverse[5][i]=
	(metric[3][i]*metric[0][i] - metric[1][i]*metric[1][i])/det[i];
    }
}


void unfused(const double metric[6][N], double inverse[6][N])
{
  double det[N];

  for(int i=0;i<N;++i)
    det[i]=metric[0][i]*metric[3][i]*metric[5][i]
      + metric[1][i]*metric[4][i]*metric[2][i]
      + metric[2][i]*metric[1][i]*metric[4][i]
      - metric[0][i]*metric[4][i]*metric[4][i]
      - metric[1][i]*metric[1][i]*metric[5][i]
      - metric[2][i]*metric[3][i]*metric[2][i];
  for(int i=0;i<N;++i)
    inverse[0][i]=
      (metric[3][i]*metric[5][i] - metric[4][i]*metric[4][i])/det[i];
  for(int i=0;i<N;++i)
    inverse[1][i]=
      (metric[2][i]*metric[4][i] - metric[1][i]*metric[5][i])/det[i];
  for(int i=0;i<N;++i)
    inverse[2][i]=
      (metric[1][i]*metric[4][i] - metric[2][i]*metric[3][i])/det[i];
  for(int i=0;i<N;++i)
    inverse[3][i]=
      (metric[0][i]*metric[5][i] - metric[2][i]*metric[2][i])/det[i];
  for(int i=0;i<N;++i)
    inverse[4][i]=
      (metric[2][i]*metric[1][i] - metric[0][i]*metric[4][i])/det[i];
  for(int i=0;i<N;++i)
    inverse[5][i]=
      (metric[3][i]*metric[0][i] - metric[1][i]*metric[1][i])/det[i];
}

#include <ctime>
#include <iostream>

using namespace std;
int main()
{
  double metric[6][N], inverse[6][N];

  for(int j=0;j<6;++j)
    for(int i=0;i<N;++i)
      metric[j][i]=1+i+j;

  const int iterations=30;

  for(int i=0;i<iterations;++i)
    {
//        fused_big(metric,inverse);
//        fused(metric,inverse);
//        unfused(metric,inverse);
    }
}

