#include<bits/stdc++.h>
// #include "stdafx.h"
#include <omp.h>
// #include<sys/time.h>

#include "matmul.h"
#include "time.h"
#include<vector>
using namespace std;

const int N = 4096;
const float seed =  0.925;


void test_omp()
{

	int A=100;
 
#pragma omp parallel for lastprivate(A)
	for(int i = 0; i<10;i++)
	{
		printf("%d\n",A);
	}
 
}


int main()
{
    
    float *a = new float[N * N];
    float *b = new float[N * N];
    float *c = new float[N * N];

    matrix_gen(a, b, N, seed);

 //   serial_matmul(a, b, c, N);
 //   parallel_matmul(a, b, c, N);
 //   splitMatMul(a, b, c, N, 32);
    splitMatMulWithSSE(a, b, c, N, 32);
  //  splitMatMulWithAVX(a, b, c, N, 32);

  // matrixMultiplyWithSSEandPartitioning(a, b, c, N, 32);

  //  vector<vector<float>> max_line(N, vector<float>(N, 0));
    float *max_line = new float[N];
    float res = 9999;

    for (int i = 0; i < N; i ++)
    {
      for (int j = 0; j < N;j ++)
      {
        max_line[i] = max(max_line[i], c[i * N + j]);
      }
    }

    for (int i = 0; i < N; i ++) res = min(res, max_line[i]);
    cout << res << endl;

    // for (int i = 0; i < N; i ++) cout << max_line[i] << ' ';

    // for (int i = 0; i < N; i ++){
    //   for (int j = 0; j < N; j ++){
    //     cout << c[i * N + j] << ' ';
    //   }
    //   cout << endl;
    // }



    return 0;
}

// grep flags -m1 /proc/cpuinfo | cut -d ":" -f 2 | tr '[:upper:]' '[:lower:]' | { read FLAGS; OPT="-march=native"; for flag in $FLAGS; do case "$flag" in "sse4_1" | "sse4_2" | "ssse3" | "fma" | "cx16" | "popcnt" | "avx" | "avx2") OPT+=" -m$flag";; esac; done; MODOPT=${OPT//_/\.}; echo "$MODOPT"; }