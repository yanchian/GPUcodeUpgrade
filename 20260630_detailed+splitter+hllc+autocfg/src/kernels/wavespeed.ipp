/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "wavespeed.hpp"

template<int blockSize>
__device__ void lastWarpReduce(volatile real* data, int tid) {
	if (blockSize >= 64) data[tid] = max(data[tid], data[tid + 32]);
	if (blockSize >= 32) data[tid] = max(data[tid], data[tid + 16]);
	if (blockSize >= 16) data[tid] = max(data[tid], data[tid + 8]);
	if (blockSize >=  8) data[tid] = max(data[tid], data[tid + 4]);
	if (blockSize >=  4) data[tid] = max(data[tid], data[tid + 2]);
	if (blockSize >=  2) data[tid] = max(data[tid], data[tid + 1]);
}

template<int blockSize>
__global__ void getFastestWaveSpeed(Mesh<GPU>::type u, Grid2D<real, GPU, 1> output) {
	int tid = threadIdx.y * blockDim.x + threadIdx.x;
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;

	__shared__ real data[512]; // maximum block size
	
	data[tid] = 0;

	if (u.active(i, j)) {
		data[tid] = getWaveSpeed(u(i, j));
	}
	__syncthreads();

	#pragma unroll
	for (int k = blockSize; k > 32; k >>= 1) {
		if (blockSize >= k * 2) {
			if (tid < k) {
				data[tid] = max(data[tid], data[tid + k]);
				__syncthreads();
			}
		}
	}
	if (tid < 32) lastWarpReduce<blockSize>(data, tid);
	
	output(blockIdx.x, blockIdx.y, 0) = data[0];
}

