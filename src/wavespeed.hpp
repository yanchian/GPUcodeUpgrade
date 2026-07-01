#pragma once

template<int blockSize>
__global__ void getFastestWaveSpeed(Mesh<GPU>::type u, Grid2D<real, GPU, 1> output);

