#include<iostream>
#include<vector>
#include<immintrin.h>
using namespace std;

float rand_float(float s){
	return 4*s*(1-s);
}


void matrix_gen(float *a,float *b,int N,float seed){
	float s=seed;
	for(int i=0;i<N*N;i++){
		s=rand_float(s);
		a[i]=s;
		s=rand_float(s);
		b[i]=s;
	}
}


// void matrix_gen(vector<vector<int>> a, vector<vector<int>> b, int N, float seed)
// {
//   //  int N = a.size();
//     float s = seed;
// 	for(int i=0;i<N;i++)
//         for (int j = 0; j < N; j ++)
//         {
//             s=rand_float(s);
//             a[i][j]=s;
//             s=rand_float(s);
//             b[i][j]=s;
// 	    }    
// }


void serial_matmul(float *a, float *b, float *c, int N)
{
    clock_t start, end;
    start = clock();
    for (int k = 0; k < N; k ++)
        for (int i = 0; i < N; i ++)
        {        
            for (int j = 0; j < N; j ++)
            {
                c[k * N + i] += a[k * N + j] * b[j * N + i];
            }      
        }
    float trace = 0;
    for (int i = 0; i < N; i ++)
        trace += c[i*N+i];
    
    cout <<"Trace:"<< trace << endl;

}



void parallel_matmul(float *a, float *b, float *c, int N)
{

    double sum = 0;
#pragma omp parallel for num_threads(16) schedule(static)
    for (int i = 0; i < N; i ++)    
        for (int j = 0; j < N; j ++)
            for (int m = 0; m < N; m ++)   
            {
            //    cout << i << ' ' << j << ' ' << m << endl;       
                c[i * N + j] += a[i * N + m] * b[m * N + j];
            }

    float trace = 0;
    for (int i = 0; i < N; i ++)
        trace += c[i*N+i];
    
    cout <<"Trace:"<< trace << endl;
}


void splitMatMul(float *a, float *b, float *c, int N, int p_size)
{
    #pragma omp parallel for num_threads(16) schedule(static)
    for (int i = 0; i < N; i += p_size )
        for (int j = 0; j < N; j += p_size )
            for (int m = 0; m < N; m += p_size)
                for (int x = i; x < i + p_size && x < N; x ++)
                    for (int y = j; y < j + p_size && y < N; y ++)
                        for (int z = m; z < m + p_size && z < N; z ++)
                            c[x * N + y] += a[x * N + z] * b[z * N + y];
    float trace = 0;
    for (int i = 0; i < N; i ++)
        trace += c[i*N+i];
    
    cout <<"Trace:"<< trace << endl;
}




// void splitMatMulWithSSE(float *a, float *b, float *c, int N, int p_size)
// {
//     __m128 pB, pA,pResult;
//     int step = 4;
   
//     for (int i = 0; i < N; i += p_size ){
//         for (int j = 0; j < N; j += p_size ){
//             for (int m = 0; m < N; m += p_size){
// #pragma omp parallel for num_threads(16) schedule(static) private(pB, pA, pResult) 
//             for (int row = m; row < m + step && row < N; row ++)
//             {
//                 for (int x = i; x < i + p_size && x < N; x += step){
//                     for (int y = j; y < j + p_size && y < N; y += step)
//                     {
//                                 pB = _mm_loadu_ps(&b[x * N + y]);
//                                 pA = _mm_set1_ps(a[row * N + x]);
//                                 pResult = _mm_mul_ps(pB, pA);

//                                 pB = _mm_loadu_ps(&b[(x + 1) * N + y]);
//                                 pA = _mm_set1_ps(a[row * N + x + 1]);
//                                 pResult = _mm_add_ps(_mm_mul_ps(pB, pA), pResult);

//                                 pB = _mm_loadu_ps(&b[(x + 2) * N + y]);
//                                 pA = _mm_set1_ps(a[row * N + x + 2]);
//                                 pResult = _mm_add_ps(_mm_mul_ps(pB, pA), pResult);

//                                 pB = _mm_loadu_ps(&b[(x + 3) * N + y]);
//                                 pA = _mm_set1_ps(a[row * N + x + 3]);
//                                 pResult = _mm_add_ps(_mm_mul_ps(pB, pA), pResult);      
//                                 _mm_storeu_ps(&c[row * N + y], _mm_add_ps(pResult, _mm_loadu_ps(&c[row * N + y])));                          
//                     }
//                 }
//             }
//             }
//         }
//     }
//     float trace = 0;
//     for (int i = 0; i < N; i ++)
//         trace += c[i*N+i];
    
//     cout <<"Trace:"<< trace << endl;
// }

void splitMatMulWithAVX(float *a, float *b, float *c, int N, int p_size)
{
  //  for (int i = 0; i < N; i ++) cout << c[i] << ' ';
    __m256 pB, pA,pResult, pC;
    int step = 8;
#pragma omp parallel for num_threads(16) schedule(static) private(pB,pA,pResult)
    for (int i = 0; i < N; i += p_size )
        for (int j = 0; j < N; j += p_size )
            for (int m = 0; m < N; m += p_size)

                for (int x = i; x < i + p_size && x < N; x += step)
                    for (int y = j; y < j + p_size && y < N; y += step)
                        for (int z = m; z < m + p_size && z < N; z += step)
                            for (int row = x; row < x + step && row < N; row ++)
                            {
                        //        cout << x << ' ' << y << ' ' << z <<  ' ' << row << endl;
                                
                                pB = _mm256_loadu_ps(&b[z * N + y]);
                                pA = _mm256_set1_ps(a[row * N + z]);
                                pResult = _mm256_mul_ps(pB, pA);

                                pB = _mm256_loadu_ps(&b[(z + 1) * N + y]);
                                pA = _mm256_set1_ps(a[row * N + z + 1]);
                                pResult = _mm256_add_ps(_mm256_mul_ps(pB, pA), pResult);

                                pB = _mm256_loadu_ps(&b[(z + 2) * N + y]);
                                pA = _mm256_set1_ps(a[row * N + z + 2]);
                                pResult = _mm256_add_ps(_mm256_mul_ps(pB, pA), pResult);

                                pB = _mm256_loadu_ps(&b[(z + 3) * N + y]);
                                pA = _mm256_set1_ps(a[row * N + z + 3]);
                                pResult = _mm256_add_ps(_mm256_mul_ps(pB, pA), pResult); 

                                pB = _mm256_loadu_ps(&b[(z + 4) * N + y]);
                                pA = _mm256_set1_ps(a[row * N + z + 4]);
                                pResult = _mm256_add_ps(_mm256_mul_ps(pB, pA), pResult);   

                                pB = _mm256_loadu_ps(&b[(z + 5) * N + y]);
                                pA = _mm256_set1_ps(a[row * N + z + 5]);
                                pResult = _mm256_add_ps(_mm256_mul_ps(pB, pA), pResult);   

                                pB = _mm256_loadu_ps(&b[(z + 6) * N + y]);
                                pA = _mm256_set1_ps(a[row * N + z + 6]);
                                pResult = _mm256_add_ps(_mm256_mul_ps(pB, pA), pResult); 

                                pB = _mm256_loadu_ps(&b[(z + 7) * N + y]);
                                pA = _mm256_set1_ps(a[row * N + z + 7]);
                                pResult = _mm256_add_ps(_mm256_mul_ps(pB, pA), pResult);       

                             //   pC = _mm256_add_ps(_mm256_loadu_ps(&c[x * N + y]), pResult);                                  
                                _mm256_storeu_ps(&c[x * N + y], _mm256_add_ps(_mm256_loadu_ps(&c[x * N + y]), pResult));                
                            }
  
    float trace = 0;
    for (int i = 0; i < N; i ++)
        trace += c[i*N+i];
    
    cout <<"Trace:"<< trace << endl;
}


void matrixMultiplyWithSSEandPartitioning(float * a, float * b, float * result, int N, int partitioning){
    __m128 pB, pA, pResult;
    int step = 4;
//#pragma omp parallel for num_threads(8) schedule(static) private(pB,pA,pResult)
    for (int BBigBlockDownMove = 0; BBigBlockDownMove < N; BBigBlockDownMove += partitioning){
        for(int BBigBlockRightMove = 0; BBigBlockRightMove < N; BBigBlockRightMove += partitioning){
            for (int row = 0; row < N; row += partitioning)
#pragma omp parallel for num_threads(16) schedule(static) private(pB,pA,pResult)
            for(int ADownMove = row; ADownMove < row + partitioning; ADownMove++){
//#pragma omp parallel for num_threads(8) schedule(static) private(pB,pA,pResult)
                    for(int BSmallBlockDownMove = BBigBlockDownMove; BSmallBlockDownMove < BBigBlockDownMove + partitioning; BSmallBlockDownMove += step){
                        for(int BSmallBlockRightMove = BBigBlockRightMove; BSmallBlockRightMove < BBigBlockRightMove + partitioning; BSmallBlockRightMove += step){
                            pB = _mm_loadu_ps(&b[BSmallBlockDownMove * N + BSmallBlockRightMove]);
                            pA = _mm_set1_ps(a[ADownMove * N + BSmallBlockDownMove]);
                            pResult = _mm_mul_ps(pB,pA);
                            pB = _mm_loadu_ps(&b[(BSmallBlockDownMove+1) * N + BSmallBlockRightMove]);
                            pA = _mm_set1_ps(a[ADownMove * N + BSmallBlockDownMove + 1]);
                            pResult = _mm_add_ps(_mm_mul_ps(pB,pA),pResult);
                            pB = _mm_loadu_ps(&b[(BSmallBlockDownMove+2) * N + BSmallBlockRightMove]);
                            pA = _mm_set1_ps(a[ADownMove * N + BSmallBlockDownMove + 2]);
                            pResult = _mm_add_ps(_mm_mul_ps(pB,pA),pResult);
                            pB = _mm_loadu_ps(&b[(BSmallBlockDownMove+3) * N + BSmallBlockRightMove]);
                            pA = _mm_set1_ps(a[ADownMove * N + BSmallBlockDownMove + 3]);
                            pResult = _mm_add_ps(_mm_mul_ps(pB,pA),pResult);
                            _mm_storeu_ps(&result[ADownMove * N + BSmallBlockRightMove], _mm_add_ps(pResult, _mm_loadu_ps(&result[ADownMove * N + BSmallBlockRightMove])));
                        }
                    }

            }
        }
    }
}