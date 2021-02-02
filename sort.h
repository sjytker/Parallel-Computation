#include<iostream>
#include<vector>
#include<immintrin.h>
using namespace std;

void sort_gen(int *d,int N,int seed){
	srand(seed);
	for(int i=0;i<N;i++){
		d[i]=rand();
	}
}

void plainRadixSort(int *d, int N) 
{
	int *buckets[10];
	int numOfBucket[10];
	int prefix = 0;
	memset(numOfBucket, 0, sizeof(numOfBucket));

//	cout << "size of numbucket : " << sizeof(numOfBucket) << endl;

#pragma omp parallel for schedule(static) num_threads(16)
	for (int i = 0; i < 10; i ++)
	{
		buckets[i] = new int[N];
		memset(buckets[i], 0, sizeof(buckets[i]));
	}

	long long divide = 1;
	int numLen = 10;
	
	for (int l = 0; l< numLen; l ++)
	{
		for (int i = 0; i < N; i ++)
		{
			
			int key = d[i] / divide % 10;
			buckets[key][numOfBucket[key]++] = d[i];
	//		cout << i << ' ' << key << endl;
		}
		divide *= 10;

        #pragma omp parallel for schedule(static) num_threads(16) private(prefix)
        for (int i = 0; i < 10; i ++)
        {
            prefix = 0;
            for (int j = 0; j < i; j ++) prefix += numOfBucket[j];
            for (int j = 0; j < numOfBucket[i]; j ++)
            {
                d[prefix + j] = buckets[i][j];
            }
        }

        for (int i = 0; i < 10; i ++) numOfBucket[i] = 0;
	}
	for (int i = 0; i < 10; i ++) delete[] buckets[i];

}


void CacheRadixSort(int *d, int N, int B, int buffer_size, int pos, int num_thread)
{

	int *thread_bucket = new int[num_thread * B]; //  count key nums
	int *c = new int[B];
	int *addr = new int[num_thread * B];
	int *buffer = new int[num_thread * B * buffer_size];
	int *buffer_cnt = new int[num_thread * B];
	int *temp = new int[N];
	memset(temp, 0, sizeof(temp));

	omp_set_num_threads(num_thread);

// different thread and bucket hash into different index
#pragma omp parallel for
	for (int i = 0; i < N; i ++)
	{
		int key = int(d[i] / pow(B, pos)) % B;
		int thread_id = omp_get_thread_num();
		thread_bucket[thread_id * B + key] ++;
	//	cout << "thread_id : " << thread_id << endl;
	}

#pragma omp barrier

	for (int j = 0; j < B; j ++)
		for (int i = 0; i < num_thread; i ++)
			c[j] += thread_bucket[i * B + j];

#pragma omp parallel for schedule(static) num_threads(num_thread)
	for (int i = 0; i < num_thread; i ++)
	{
		for (int j = 0; j < B; j ++)
		{
			for (int k = 0; k < j; k ++) addr[i * B + j] += c[k];
			for (int k = 0; k < i; k ++) addr[i * B + j] += thread_bucket[k * B + j];
		}
	}

	int cnt = 0;
	memset(buffer_cnt, 0, sizeof(buffer_cnt));
#pragma omp parallel for
	for (int i = 0; i < N; i ++)
	{
		int key = int(d[i] / pow(B, pos)) % B;
		int thread_id = omp_get_thread_num();
		buffer[ thread_id * B + key + buffer_cnt[thread_id * B + key]++ ] = d[i];

		if (buffer_cnt[thread_id * B + key] == buffer_size)
		{
			for (int j = 0; j < buffer_size; j ++)
				temp[addr[thread_id * B + key] ++] = buffer[thread_id * B + key + j];
	//		addr[thread_id * B + key] += buffer_size;
			buffer_cnt[thread_id * B + key] = 0;
		}
	}

// #pragma omp parallel for schedule(static) num_threads(num_thread)
	for (int k = 0; k < num_thread; k ++)
	{
		int thread_id = omp_get_thread_num();
		for (int j = 0;j < B; j ++)
		{
			int t = buffer_cnt[k * B + j];
			if (t)
			{
				while (t--) temp[addr[k * B + j] ++] = buffer[k * B + j];
			}
		}
	}	

	memcpy(d, temp, sizeof(d));
	return ;
}



void creat_hash(int* s, int N, int B, int* D, int pos, int thread_num,int l) {

	
	int** p;
	int* c = new int[B]();
	int d, h;
	int** b_thread;//每锟斤拷锟竭筹拷锟斤拷锟斤拷锟斤拷锟紹锟斤拷桶锟叫碉拷元锟斤拷锟斤拷
	b_thread = new int* [thread_num];
	for (int i = 0; i < thread_num; i++) b_thread[i] = new int[B]();
	p = new int* [thread_num];
	for (int i = 0; i < thread_num; i++) p[i] = new int[B]();


	omp_set_num_threads(thread_num);
int thread_id;
#pragma omp parallel for private(thread_id,d,h) 
	for (int i = 0; i < N; i++) {
		 thread_id = omp_get_thread_num();
		d = s[i] / pow(B, pos);
		h = d % B;

		b_thread[thread_id][h]++;

	}

	for (int j = 0; j < B; j++) {
		for (int i = 0; i < thread_num; i++) {
			c[j] += b_thread[i][j];
		}
	}

	for (int i = 1; i < thread_num; i++) {
		for (int j = 1; j < B; j++) {
		//	for (int k = 0; k < j; k++) p[i][j] += c[k];
		//	for (int k = 0; k < i; k++) p[i][j] += b_thread[k][j];
            p[i][j] += p[i][j - 1] + c[j - 1];
            p[i][j] += p[i - 1][j] + b_thread[i - 1][j]
		}
	}

	vector<vector<vector<int>>> b(thread_num, vector<vector<int> >(B));//锟捷达拷锟斤拷锟斤拷锟矫匡拷锟酵帮拷锟斤拷锟�
#pragma omp parallel for private(thread_id,d,h) 
	for (int i = 0; i < N; i++) {
		int thread_id = omp_get_thread_num();
		d = s[i] / pow(B, pos);
		h = d % B;
		b[thread_id][h].push_back(s[i]);
		if (b[thread_id][h].size() == l) {
			for (int j = 0; j < l; j++) {
				D[p[thread_id][h]+j] = b[thread_id][h][j];
			}
			p[thread_id][h] = p[thread_id][h] + l;
			b[thread_id][h].clear();
		}

		
		//c[h]++;
		//cout << "thread_id:" << thread_id << "  i:" << i << endl;

	}

	for(int k=0;k<thread_num;k++){
		for (int i = 0; i < B; i++) {
			if (!b[k][i].empty()) {
				int len = b[k][i].size();
				for (int j = 0; j < len; j++) {
					D[p[k][i]+j] = b[k][i][j];
				}
				p[k][i] = p[k][i] + len;
			}
		}
	}


	delete[]c;

	for (int i = 0; i < thread_num; i++)
		delete[] b_thread[i];
	delete[] b_thread;
	for (int i = 0; i < thread_num; i++)
		delete[] p[i];
	delete[] p;
}


void merge_sort(int *d, int left, int right, int N)
{
	if (left == right) return;
	int mid = (left + right) / 2;
	merge_sort(d, left, mid, N);
	merge_sort(d, mid + 1, right, N);
	int *tmp = new int[N];
	
	int k = left;
	int l = left, r = mid + 1;
	for (; l <= mid && r <= right;)
	{
		if (d[l] <= d[r]) tmp[k++] = d[l++];
		else tmp[k++] = d[r++];
	}
	while (l <= mid) tmp[k++] = d[l++];
	while (r <= right) tmp[k++] = d[r++];
	for (k = left; k <= right;k++) d[k] = tmp[k];
	delete[]tmp;
}

void mergeSort(int * array, int left, int right, int N){
    if(left == right)return;
    int mid = (left + right)/2;
    mergeSort(array, left, mid, N);
    mergeSort(array, mid+1, right, N);
    int *tmpArray = new int [N];
    memcpy(tmpArray, array, N * sizeof(int));
    int i1 = left, i2 = mid+1;
    for(int curr = left; curr<=right; curr++){
        if(i1 == mid + 1)
            array[curr] = tmpArray[i2++];
        else if(i2>right)
            array[curr] = tmpArray[i1++];
        else if(tmpArray[i1] < tmpArray[i2])
            array[curr] = tmpArray[i1++];
        else array[curr] = tmpArray[i2++];
    }
    delete[]tmpArray;
}

void parallel_mergeSort(int *d, int N)
{
	int step = 2;
	int *tmp = new int[N];
	int left, right, mid, i1, i2, k;

	cout << " in pSort " << endl;
	while (step <= N)
	{
#pragma omp parallel for num_threads(16) private(left, right, mid, i1, i2, k)
		for (int i = 0; i <= N; i += step)
		{
			left = i;
			right = i + step - 1;
			mid = right + (right - left) >> 1;

			cout << "i left right : " << i << ' ' << left << ' ' << right << ' ' << step << endl;

			i1 = left, i2 = mid + 1;
			k = left;
			for (; i1 <= mid && i2 <= right;)
			{
				if (d[i1] <= d[i2]) tmp[k++] = d[i1++];
				else tmp[k++] = d[i2++];
			}
			while (i1 <= mid) tmp[k++] = d[i1++];
			while (i2 <= right) tmp[k++] = d[i2++];

			for (k = left; k <= right; k ++) d[k] = tmp[k];
		}
		step *= 2;
		cout << "step N : " << step << ' ' << N << endl;
	}
	cout << "finish" << endl;
//	delete[] tmp;
}

void mergeSortParallel(int *array, int N){
    unsigned int step = 2;
    int maximun = 2147483647;
    int  right, mid, i1, i2;
    int *tmpArray_1, *tmpArray_2;
    int bit = N;
    double log2Result = log2(N);
    double ceilingResult = ceil(log2Result);

    if(log2Result != ceilingResult){
        bit = pow(2, ceilingResult);
    }

    tmpArray_1 = new int[bit];

    if(log2Result != ceilingResult){
        tmpArray_2 = new int[bit];
#pragma omp parallel for num_threads(16)
        for(int i = 0; i < bit; i++){
            tmpArray_1[i] = maximun;
        }
        memcpy(tmpArray_1, array, N * sizeof(int));
        memcpy(tmpArray_2, tmpArray_1, bit * sizeof(int));
        while(step < bit)
        {

            memcpy(tmpArray_1, tmpArray_2, bit * sizeof(int));
#pragma omp parallel for num_threads(16) private( right, mid, i1, i2)
            for(unsigned int i = 0; i < bit; i += step){
                int left = i;
                right = i + step -1;
                mid = (left + right) / 2;
                i1 = left;
                i2 = mid + 1;
                for(unsigned int curr = left; curr<=right; curr++){
                    if(i1 == mid + 1)
                        tmpArray_2[curr] = tmpArray_1[i2++];
                    else if(i2>right)
                        tmpArray_2[curr] = tmpArray_1[i1++];
                    else if(tmpArray_1[i1] < tmpArray_1[i2])
                        tmpArray_2[curr] = tmpArray_1[i1++];
                    else tmpArray_2[curr] = tmpArray_1[i2++];
                }
            }
            step *= 2;
        }

        memcpy(array, tmpArray_2, N * sizeof(int));
        delete[] tmpArray_2;
    }
    else{
        while(step <= N)
        {
            memcpy(tmpArray_1, array, N * sizeof(int));
#pragma omp parallel for private( right, mid, i1, i2)
            for(unsigned int i = 0; i < N; i += step){
                int left = i;
                right = i + step -1;
                mid = (left + right) / 2;
                i1 = left;
                i2 = mid + 1;
                for(unsigned int curr = left; curr<=right; curr++){
                    if(i1 == mid + 1)
                        array[curr] = tmpArray_1[i2++];
                    else if(i2>right)
                        array[curr] = tmpArray_1[i1++];
                    else if(tmpArray_1[i1] < tmpArray_1[i2])
                        array[curr] = tmpArray_1[i1++];
                    else array[curr] = tmpArray_1[i2++];
                }
            }
            step *= 2;
        }
    }
    delete[] tmpArray_1;
}

void mergeForSSE(int *tmpArray_1, int * tmpArray_2,int bit, int step){
    int left, right, mid, i1, i2;
    if(step == 2){
#pragma omp parallel for num_threads(16) private(left, right, mid, i1, i2)
        for(unsigned int i = 0; i < bit; i += step){
            left = i;
            right = i + step -1;
            mid = (left + right) / 2;
            i1 = left;
            i2 = mid + 1;
            for(unsigned int curr = left; curr<=right; curr++){
                if(i1 == mid + 1)
                    tmpArray_2[curr] = tmpArray_1[i2++];
                else if(i2>right)
                    tmpArray_2[curr] = tmpArray_1[i1++];
                else if(tmpArray_1[i1] < tmpArray_1[i2])
                    tmpArray_2[curr] = tmpArray_1[i1++];
                else tmpArray_2[curr] = tmpArray_1[i2++];
            }
        }
    }
    else{
#pragma omp parallel for num_threads(16) private(left, right, mid, i1, i2)
        for(unsigned int i = 0; i < bit; i += step){
            left = i;
            right = i + step -1;
            mid = (left + right) / 2;
            i1 = left;
            i2 = mid + 1;
            if(!(i / step % 2)){
                for(unsigned int curr = left; curr<=right; curr++){
                    if(i1 == mid + 1)
                        tmpArray_2[curr] = tmpArray_1[i2++];
                    else if(i2>right)
                        tmpArray_2[curr] = tmpArray_1[i1++];
                    else if(tmpArray_1[i1] < tmpArray_1[i2])
                        tmpArray_2[curr] = tmpArray_1[i1++];
                    else tmpArray_2[curr] = tmpArray_1[i2++];
                }
            }else{
                for(unsigned int curr = left; curr<=right; curr++){
                    if(i1 == mid + 1)
                        tmpArray_2[right - curr + left] = tmpArray_1[i2++];
                    else if(i2>right)
                        tmpArray_2[right - curr + left] = tmpArray_1[i1++];
                    else if(tmpArray_1[i1] < tmpArray_1[i2])
                        tmpArray_2[right - curr + left] = tmpArray_1[i1++];
                    else tmpArray_2[right - curr + left] = tmpArray_1[i2++];
                }
            }

        }
    }
}

void bitonicMergeNetworkForSSE(int *tmpArray_1, int * tmpArray_2,int bit, int step){

    int left, right, mid, end;
    __m128i srciA, srciB, L_1, H_1, dstiA, dstiB;
    __m128i midL_1, midH_1;
    int one = 0xFFFFFFFF, zero = 0X00000000;
//    __m128i loadIndex = _mm_set_epi32(0, 1, 2, 3);
    __m128i storeMask = _mm_set_epi32(one, one, one, one);
    __m128i mask13ZERO = _mm_set_epi32(one, zero, one, zero);
    __m128i mask24ZERO = _mm_set_epi32(zero, one, zero, one);
    int tmp;
//#pragma omp parallel for private(srciA, srciB, L_1, H_1, dstiA, dstiB, midL_1, midH_1, left, right, mid, end) firstprivate(step, bit)
    for(unsigned int i = 0; i < bit; i += step) {
        left = i;
        end = i + step - 1;
        mid = (left + end) / 2;
        right = mid + 1;
        for(int k = 0; k < log2(step); k++){
#pragma omp parallel for num_threads(16) private(srciA, srciB, L_1, H_1, dstiA, dstiB, midL_1, midH_1) firstprivate(step, bit)
            for(unsigned int j = 0; j < step / 2; j += 4) {

                srciA = _mm_set_epi32(tmpArray_1[left + j], tmpArray_1[left + j + 1], tmpArray_1[left + j + 2], tmpArray_1[left + j + 3]);
                srciB = _mm_set_epi32(tmpArray_1[right + j], tmpArray_1[right + j + 1], tmpArray_1[right + j + 2], tmpArray_1[right + j + 3]);
                H_1 = _mm_max_epi32(srciA, srciB);
                L_1 = _mm_min_epi32(srciA, srciB);
                midL_1 = _mm_shuffle_epi32(L_1, 0x27);
                midH_1 = _mm_shuffle_epi32(H_1, 0x8D);
                midL_1 = _mm_and_si128(midL_1, mask24ZERO);
                midH_1 = _mm_and_si128(midH_1, mask13ZERO);
                dstiA = _mm_add_epi32(midH_1, midL_1);

                midL_1 = _mm_shuffle_epi32(L_1, 0x8D);
                midH_1 = _mm_shuffle_epi32(H_1, 0x27);
                midL_1 = _mm_and_si128(midL_1, mask24ZERO);
                midH_1 = _mm_and_si128(midH_1, mask13ZERO);
                dstiB = _mm_add_epi32(midH_1, midL_1);

                _mm_maskstore_epi32(&tmpArray_2[left + j * 2], storeMask, dstiA);
                _mm_maskstore_epi32(&tmpArray_2[left + j * 2 + 4], storeMask, dstiB);


            }
            memcpy(&tmpArray_1[left], &tmpArray_2[left], step * sizeof(int));
        }



        if (i / step % 2) {
//#pragma omp parallel for private(tmp)
            for(int m = 0; m < step / 2; m++){
                tmp = tmpArray_2[left + m];
                tmpArray_2[left + m] = tmpArray_2[end - m];
                tmpArray_2[end - m] = tmp;
            }
        }
    }
}

void bitonicMergeNetworkForSSEv2(int *tmpArray_1, int * tmpArray_2,int bit, int step){

    int left, right, mid, end;
    __m128i srciA, srciB, L_1, H_1, dstiA, dstiB;
    __m128i midL_1, midH_1;
    int one = 0xFFFFFFFF, zero = 0X00000000;
//    __m128i loadIndex = _mm_set_epi32(0, 1, 2, 3);
    __m128i storeMask = _mm_set_epi32(one, one, one, one);
    __m128i mask13ZERO = _mm_set_epi32(one, zero, one, zero);
    __m128i mask24ZERO = _mm_set_epi32(zero, one, zero, one);
    int tmp;
    for(int k = 0;k < log2(step); k++){
#pragma omp parallel for num_threads(16) private(srciA, srciB, L_1, H_1, dstiA, dstiB, midL_1, midH_1, left, right, mid, end) firstprivate(step, bit)
        for(unsigned int i = 0; i < bit; i += step) {
            left = i;
            end = i + step - 1;
            mid = (left + end) / 2;
            right = mid + 1;
//#pragma omp parallel for private(srciA, srciB, L_1, H_1, dstiA, dstiB, midL_1, midH_1) firstprivate(step, bit)
            for(unsigned int j = 0; j < step / 2; j += 4) {

                srciA = _mm_set_epi32(tmpArray_1[left + j], tmpArray_1[left + j + 1], tmpArray_1[left + j + 2], tmpArray_1[left + j + 3]);
                srciB = _mm_set_epi32(tmpArray_1[right + j], tmpArray_1[right + j + 1], tmpArray_1[right + j + 2], tmpArray_1[right + j + 3]);
                H_1 = _mm_max_epi32(srciA, srciB);
                L_1 = _mm_min_epi32(srciA, srciB);
                midL_1 = _mm_shuffle_epi32(L_1, 0x27);
                midH_1 = _mm_shuffle_epi32(H_1, 0x8D);
                midL_1 = _mm_and_si128(midL_1, mask24ZERO);
                midH_1 = _mm_and_si128(midH_1, mask13ZERO);
                dstiA = _mm_add_epi32(midH_1, midL_1);

                midL_1 = _mm_shuffle_epi32(L_1, 0x8D);
                midH_1 = _mm_shuffle_epi32(H_1, 0x27);
                midL_1 = _mm_and_si128(midL_1, mask24ZERO);
                midH_1 = _mm_and_si128(midH_1, mask13ZERO);
                dstiB = _mm_add_epi32(midH_1, midL_1);

                _mm_maskstore_epi32(&tmpArray_2[left + j * 2], storeMask, dstiA);
                _mm_maskstore_epi32(&tmpArray_2[left + j * 2 + 4], storeMask, dstiB);
            }
        }
        memcpy(tmpArray_1, tmpArray_2, bit * sizeof(int));
    }
    for(unsigned int i = 0; i < bit; i += step) {
        left = i;
        end = i + step - 1;
        if (i / step % 2) {
            for(int m = 0; m < step / 2; m++){
                tmp = tmpArray_2[left + m];
                tmpArray_2[left + m] = tmpArray_2[end - m];
                tmpArray_2[end - m] = tmp;
            }
        }
    }
}

void bitonicSortBaseOnSSE(int * array, int N){
    unsigned int step = 2;
    int *tmpArray_1 = new int[N];

    memcpy(tmpArray_1, array, N * sizeof(int));
    while(step <= N)
    {

        memcpy(tmpArray_1, array, N * sizeof(int));

        if(step == 2 || step == 4){
            mergeForSSE(tmpArray_1, array, N, step);
        }
        else{
            bitonicMergeNetworkForSSEv2(tmpArray_1, array, N, step);
        }
        step *= 2;
    }

    delete[] tmpArray_1;
}
