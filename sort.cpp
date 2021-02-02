#include<bits/stdc++.h>
#include <omp.h>
#include<sys/time.h>
#include "sort.h"
using namespace std;

const int N = 1000000000;
int seed = 99;

void calcLen(int x, int &l, int &n)
{
    int res = 0;
    while (x)
    {
        res ++;
        n = x % 10;
        x /= 10;
    }
    l = res;
}

int main()
{
 //   int* d = (int*)malloc(N * sizeof(int));
    int *d = new int[N];

    timeval t_start, tp_start, t_end;
    gettimeofday(&t_start, NULL);

    sort_gen(d, N, seed);

    const int K = 32;
    const int D0 = 257;
    const int thread = 8;

    // for (int K = 16; K <= 64; K *= 2)
    //     for (int D0 = 255; D0 <= 257; D0 ++)

    gettimeofday(&tp_start, NULL);

  // sort(d, d + N);   //367s, 1075275326
 
 //   plainRadixSort(d, N);    // 202s
    // creat_hash(d1, N, D0, d2, 0, thread, K);
    // creat_hash(d2, N, D0, d1, 1, thread, K);
    // creat_hash(d1, N, D0, d2, 2, thread, K);
    // creat_hash(d2, N, D0, d1, 3, thread, K);

 //    merge_sort(d, 0, N - 1, N);

 //   parallel_mergeSort(d, N);

 //   mergeSortParallel(d, N);

    bitonicSortBaseOnSSE(d, N);

    cout << d[N / 2] << endl;

    gettimeofday( &t_end, NULL);
    double delta_t = (t_end.tv_sec-t_start.tv_sec) + 
                    (t_end.tv_usec-t_start.tv_usec) / 1000000.0;

    double delta_tp = (t_end.tv_sec-tp_start.tv_sec) + 
                    (t_end.tv_usec-tp_start.tv_usec) / 1000000.0;

    // cout << "sort time : " << delta_tp  << "s" << ' ' << K << ' ' << D0 << endl;
    cout << "sort time : " << delta_tp << "s" << endl;
    cout << "all time : " << delta_t  << "s" << endl;


    return 0;
}