#pragma once

template<int blockDimX, int blockDimY, bool X, bool SET>
__global__ void getDiffusionFluxesKernel(Mesh<GPU>::type u, Mesh<GPU>::type flux, const real dt);


