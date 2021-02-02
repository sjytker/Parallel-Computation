

研一并行计算课程实验



#### 1. 并行矩阵乘法

- 设置矩阵分块以保证 cache 命中
- simd 指令并行化，包括 sse, avx
- openmp fork-join





#### 2. 并行基数排序

- 根据 L1 cache 64 KB 设置合适的 buffer = 32 长度数组，基数大小 B = 256。这样 32 x 256 x 4 x 2 = 64 KB，左右数组值恰好利用完毕 L1 cache。注意， L1 cache 包括 icache, dcache，也就是存指令和数据的，在 windows 下直接查看会显示 128 KB，但在 linux 命令行中查询会分开显示 64 KB

- fork-join 并行哈希当前基数位，统计每个线程 i 每一位 j 有多少个数 `b_thread[thread_id][key]` 

- 对 `b_thread[thread_id][key]` 计算线程 i，基数 j 的前缀和数组 `p[i][j]`，这用于给下一步并行存入目标数组提供地址。

  ```c
  	for (int i = 1; i < thread_num; i++) {
  		for (int j = 1; j < B; j++) {
              p[i][j] += p[i][j - 1] + c[j - 1];
              p[i][j] += p[i - 1][j] + b_thread[i - 1][j]
  		}
  	}
  ```

- 最后要先把原数组并行哈希到一个临时数组，当临时数组长度到达预设（这里是32），再把它写入目标数组。这是因为临时数组和目标数组存储在不同的地方，如果每次把一个哈希结果逐个写入目标数组，cache miss 会非常严重

- 完毕后再把剩余数（未达到 buffer size 的）存入目标数组



#### 3. 并行归并排序

自底向上归并即可，天然地划分出了 k 个区间，直接并行化