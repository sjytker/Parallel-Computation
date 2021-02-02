#include<bits/stdc++.h>
#include <omp.h>
#include<sys/time.h>
#include "matmul.h"
#include<vector>
using namespace std;

const int N = 8192;
// float a[N][N], b[N][N], c[N][N];


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
  //  test_omp();
    
    float seed = 0.3;
    float *a = new float[N * N];
    float *b = new float[N * N];
    float *c = new float[N * N];
    timeval t_start, tp_start, t_end;
    gettimeofday(&t_start, NULL);
    // vector<vector<int>> a(N, vector<int>(N, 0));
    // vector<vector<int>> b(N, vector<int>(N, 0));
    // vector<vector<int>> c(N, vector<int>(N, 0));

    matrix_gen(a, b, N, seed);

    gettimeofday(&tp_start, NULL);
  //  serial_matmul(a, b, c, N);
  //  parallel_matmul(a, b, c, N);
  //  splitMatMul(a, b, c, N, 32);
  //  splitMatMulWithSSE(a, b, c, N, 32);
    splitMatMulWithAVX(a, b, c, N, 32);
 //   matrixMultiplyWithSSE(a, b, c, N);

    gettimeofday( &t_end, NULL);
    double delta_t = (t_end.tv_sec-t_start.tv_sec) + 
                    (t_end.tv_usec-t_start.tv_usec) / 1000000.0;

    double delta_tp = (t_end.tv_sec-tp_start.tv_sec) + 
                    (t_end.tv_usec-tp_start.tv_usec) / 1000000.0;

    cout << "parrallel time : " << delta_tp  << "s" << endl;
    cout << "elapsed time : " << delta_t  << "s" << endl;
    return 0;
}

// grep flags -m1 /proc/cpuinfo | cut -d ":" -f 2 | tr '[:upper:]' '[:lower:]' | { read FLAGS; OPT="-march=native"; for flag in $FLAGS; do case "$flag" in "sse4_1" | "sse4_2" | "ssse3" | "fma" | "cx16" | "popcnt" | "avx" | "avx2") OPT+=" -m$flag";; esac; done; MODOPT=${OPT//_/\.}; echo "$MODOPT"; }