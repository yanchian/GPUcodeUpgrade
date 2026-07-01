/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/

__global__ void pMaxKernel(Mesh<GPU>::type u, Grid2D<real, GPU, 1> p) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;

	int offset = (u.time() * 2845) / u.dx();

	if (u.active(i, j) && i + offset < p.Nx()) {
		if (u(i, j, PRESSURE) > p(i + offset, j, 0)) {
			p(i + offset, j, 0) = u(i, j, PRESSURE);
		}
	}
}

int pMax(Mesh<GPU>::type u, Grid2D<real, GPU, 1> p) {
	dim3 blockDim(16, 16);
	dim3 gridDim = u.activeCover(blockDim);

	pMaxKernel<<<gridDim, blockDim>>>(u, p);

	return (u.time() * 2845) / u.dx() > p.Nx();
}

