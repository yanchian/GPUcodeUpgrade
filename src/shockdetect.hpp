#pragma once

__global__ void shockDetectKernel(Mesh<GPU>::type u);

__global__ void shockBurnKernel(Mesh<GPU>::type u);

