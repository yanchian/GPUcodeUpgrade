__global__ void shockDetectKernel(Mesh<GPU>::type u) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;

	const real pThreshold = 2;

	if (u.active(i, j)) {
		const real pLeft = u(i - 2, j, PRESSURE);
		const real pRight = u(i + 2, j, PRESSURE);
		const real pUp = u(i, j - 2, PRESSURE);
		const real pDown = u(i, j + 2, PRESSURE);

		if (abs(pLeft - pRight) / min(pLeft, pRight) > pThreshold || abs(pUp - pDown) / min(pUp, pDown) > pThreshold) {
			u(i, j, ISSHOCK) = 1.0;
		} else {
			u(i, j, ISSHOCK) = 0.0;
		}
	}
}

__global__ void shockBurnKernel(Mesh<GPU>::type u) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;

	const real pThreshold = 2;

	if (u.active(i, j)) {
		if (u(i, j, ISSHOCK)) {
			for (int ii = -1; ii <= 1; ii++) {
				for (int jj = -1; jj <= 1; jj++) {
					u(i + ii, j + jj, ISSHOCK) = 1.0;
				}
			}
		}
	}
}

