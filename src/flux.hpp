#pragma once

__device__ __host__ void conservativeToPrimitive(const Cell u, Cell p);

__device__ __host__ void conservativeToPrimitiveInPlace(Cell u);

__device__ __host__ void primitiveToConservativeInPlace(Cell p);

__device__ __host__ real soundSpeed(const Cell u);

__device__ __host__ real soundSpeedPrimitive(const Cell p);

__device__ real getWaveSpeed(const Cell& u);

template<bool X>
__device__ void conservativeToFlux(const Cell u, Cell f);

__device__ void extrapolate(Cell uL, const Cell uC, Cell uR);

template<int X>
__device__ void primitiveToFlux(const Cell p, Cell f);

template<int blockDimX, int blockDimY, bool X, bool SET>
__global__ void getMUSCLFluxes(Mesh<GPU>::type u, Mesh<GPU>::type flux, const real dt);

template<bool X, bool SET>
__global__ void addFluxesKernel(Mesh<GPU>::type u, Mesh<GPU>::type flux);

__global__ void addSemiFluxesKernel(Mesh<GPU>::type u, Mesh<GPU>::type flux);

