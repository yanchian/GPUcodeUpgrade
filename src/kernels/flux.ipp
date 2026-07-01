/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "flux.hpp"
//#include "kernels/HLLE.ipp"
//#include "kernels/HLLC.ipp"

__device__ __host__ void conservativeToPrimitive(const Cell u, Cell p) {
  p[DENSITY]   = u[DENSITY];
  p[XVELOCITY] = u[XMOMENTUM] / u[DENSITY];
  p[YVELOCITY] = u[YMOMENTUM] / u[DENSITY];
  p[PRESSURE]  = (u[ENERGY] - 0.5 * p[DENSITY] * (p[XVELOCITY] * p[XVELOCITY] + p[YVELOCITY] * p[YVELOCITY])) * (gamma(u) - 1.0) - gamma(u) * p0(u);

  // Species mass fractions are stored directly (no density multiplication)
  for (int k = 0; k < NUM_SPECIES - 1; k++) {
    p[SPECIES_INDEX(k)] = u[SPECIES_INDEX(k)];
  }

  for (int i = CONSERVATIVE_VARIABLES; i < NUMBER_VARIABLES; i++) {
    p[i] = u[i];
  }
}

__device__ __host__ void conservativeToPrimitiveInPlace(Cell u) {
  u[XVELOCITY] = u[XMOMENTUM] / u[DENSITY];
  u[YVELOCITY] = u[YMOMENTUM] / u[DENSITY];
  //u.momentum() = u.momentum() / u.density();
  //u[PRESSURE]  = (u[ENERGY] - 0.5 * u.density() * norm2(u.momentum())) * (gamma(u) - 1.0) - gamma(u) * p0(u);
  u[PRESSURE]  = (u[ENERGY] - 0.5 * u[DENSITY] * (u[XVELOCITY] * u[XVELOCITY] + u[YVELOCITY] * u[YVELOCITY])) * (gamma(u) - 1.0) - gamma(u) * p0(u);
  // Species are already stored as mass fractions (no conversion needed)
}

__device__ __host__ void primitiveToConservative(const Cell p, Cell u) {
  u[DENSITY]   = p[DENSITY];
  u[XMOMENTUM] = p[XVELOCITY] * p[DENSITY];
  u[YMOMENTUM] = p[YVELOCITY] * p[DENSITY];
  u[ENERGY]    = (p[PRESSURE] + gamma(p) * p0(p)) / (gamma(p) - 1.0) + 0.5 * (p[XVELOCITY] * p[XVELOCITY] + p[YVELOCITY] * p[YVELOCITY]) * p[DENSITY];
  for (int k = 0; k < NUM_SPECIES - 1; k++) {
    u[SPECIES_INDEX(k)] = p[SPECIES_INDEX(k)];
  }
  for (int i = CONSERVATIVE_VARIABLES; i < NUMBER_VARIABLES; i++) {
    u[i] = p[i];
  }
}

__device__ __host__ void primitiveToConservativeInPlace(Cell p) {
  p[XMOMENTUM] = p[XVELOCITY] * p[DENSITY];
  p[YMOMENTUM] = p[YVELOCITY] * p[DENSITY];
  p[ENERGY]    = (p[PRESSURE] + gamma(p) * p0(p)) / (gamma(p) - 1.0) + 0.5 * (p[XMOMENTUM] * p[XMOMENTUM] + p[YMOMENTUM] * p[YMOMENTUM]) / p[DENSITY];
  // Species already stored as mass fractions (no conversion needed)
}

__device__ __host__ real soundSpeed(const Cell u) {
  const real pressure = (u[ENERGY] - 0.5 * (u[XMOMENTUM] * u[XMOMENTUM] + u[YMOMENTUM] * u[YMOMENTUM]) / u[DENSITY]) * (gamma(u) - 1.0) - gamma(u) * p0(u);
  real c_3 = gamma(u) * (pressure + p0(u)) / u[DENSITY];
  if (c_3 < 1e-9){
	  c_3 = 1.0;
  }
  return sqrt(c_3);
//  return sqrt(gamma(u) * (pressure + p0(u)) / u[DENSITY]);
}

__device__ __host__ real soundSpeedPrimitive(const Cell p) {
	real c_2 = gamma(p) * (p[PRESSURE] + p0(p)) / p[DENSITY];
	if ( c_2 < 1e-9) {
		c_2 = 1.0;
	}
  return sqrt(c_2);
  // return sqrt(gamma(p) * (p[PRESSURE] + p0(p)) / p[DENSITY]);
}
__device__ real getWaveSpeed(const Cell& u) {
#ifdef GHOST
  if (u[PHI] >= 0) {
    return 0;
  } else
#endif
{
  real p[NUMBER_VARIABLES];
    conservativeToPrimitive(u, p);
	const real c = soundSpeedPrimitive(p);
	return c + sqrt(p[XVELOCITY] * p[XVELOCITY] + p[YVELOCITY] * p[YVELOCITY]);
}
}

template<bool X>
__device__ void conservativeToFlux(const Cell u, Cell f) {
  const real pressure = (u[ENERGY] - (u[XMOMENTUM] * u[XMOMENTUM] + u[YMOMENTUM] * u[YMOMENTUM]) / (2.0 * u[DENSITY])) * (gamma(u) - 1) - gamma(u) * p0(u);

  if (X) {
    f[DENSITYFLUX]   = u[XMOMENTUM];
    f[XMOMENTUMFLUX] = u[XMOMENTUM] * u[XMOMENTUM] / u[DENSITY] + pressure;
    f[YMOMENTUMFLUX] = u[XMOMENTUM] * u[YMOMENTUM] / u[DENSITY];
    f[ENERGYFLUX]    = u[XMOMENTUM] / u[DENSITY] * (u[ENERGY] + pressure);
    for (int k = 0; k < NUM_SPECIES - 1; k++) {
      f[SPECIES_INDEX(k)] = u[SPECIES_INDEX(k)] * u[XMOMENTUM] / u[DENSITY];
    }
    for (int i = CONSERVATIVE_VARIABLES; i < CONSERVATIVE_VARIABLES + NONCONSERVATIVE_VARIABLES; i++) {
      f[i]         = u[XMOMENTUM] / u[DENSITY] * u[i];
    }
  } else {
    f[DENSITYFLUX]   = u[YMOMENTUM];
    f[XMOMENTUMFLUX] = u[XMOMENTUM] * u[YMOMENTUM] / u[DENSITY];
    f[YMOMENTUMFLUX] = u[YMOMENTUM] * u[YMOMENTUM] / u[DENSITY] + pressure;
    f[ENERGYFLUX]    = u[YMOMENTUM] / u[DENSITY] * (u[ENERGY] + pressure);
    for (int k = 0; k < NUM_SPECIES - 1; k++) {
      f[SPECIES_INDEX(k)] = u[SPECIES_INDEX(k)] * u[YMOMENTUM] / u[DENSITY];
    }
    for (int i = CONSERVATIVE_VARIABLES; i < CONSERVATIVE_VARIABLES + NONCONSERVATIVE_VARIABLES; i++) {
      f[i]         = u[YMOMENTUM] / u[DENSITY] * u[i];
    }
  }
}

__device__ real limiter(real r) {
  // van Leer
  if (r < 0) return 0;
  return 2.0 * min(1.0, r) / (1.0 + r);
}

__device__ void extrapolate(Cell uL, const Cell uC, Cell uR) {
  for (int k = 0; k < CONSERVATIVE_VARIABLES + NONCONSERVATIVE_VARIABLES; k++) {
    // compute \Delta_i, \Delta{i+1}
    real delta_minus  = uC[k] - uL[k];
    real delta_plus   = uR[k] - uC[k];

    // cap them from below to some small constant
    delta_minus = copysign(max(1e-8, abs(delta_minus)), delta_minus);
    delta_plus  = copysign(max(1e-8, abs(delta_plus)), delta_plus);

    const real xi = limiter(delta_minus / delta_plus);

    // now limit the central difference approximation to \Delta (c = 0)
    const real gradient = xi * 0.25 * (uR[k] - uL[k]);

    // update extrapolations
    uL[k] = uC[k] - gradient;
    uR[k] = uC[k] + gradient;
  }
}

template<int X>
__device__ void primitiveToFlux(const Cell p, Cell f) {
  if (X) {
    f[DENSITYFLUX]   = p[DENSITY] * p[XVELOCITY];
    f[XMOMENTUMFLUX] = p[DENSITY] * p[XVELOCITY] * p[XVELOCITY] + p[PRESSURE];
    f[YMOMENTUMFLUX] = p[DENSITY] * p[XVELOCITY] * p[YVELOCITY];
    f[ENERGYFLUX]    = p[XVELOCITY] * ((p[PRESSURE] + p0(p)) * gamma(p) / (gamma(p) - 1.0) + 0.5 * p[DENSITY] * (p[XVELOCITY] * p[XVELOCITY] + p[YVELOCITY] * p[YVELOCITY]));
    for (int k = 0; k < NUM_SPECIES - 1; k++) {
      f[SPECIES_INDEX(k)] = p[SPECIES_INDEX(k)] * p[XVELOCITY];
    }
    for (int i = CONSERVATIVE_VARIABLES; i < CONSERVATIVE_VARIABLES + NONCONSERVATIVE_VARIABLES; i++) {
      f[i]         = p[XVELOCITY] * p[i];
    }
  } else {
    f[DENSITYFLUX]   = p[DENSITY] * p[YVELOCITY];
    f[XMOMENTUMFLUX] = p[DENSITY] * p[XVELOCITY] * p[YVELOCITY];
    f[YMOMENTUMFLUX] = p[DENSITY] * p[YVELOCITY] * p[YVELOCITY] + p[PRESSURE];
    f[ENERGYFLUX]    = p[YVELOCITY] * ((p[PRESSURE] + p0(p)) * gamma(p) / (gamma(p) - 1.0) + 0.5 * p[DENSITY] * (p[XVELOCITY] * p[XVELOCITY] + p[YVELOCITY] * p[YVELOCITY]));
    for (int k = 0; k < NUM_SPECIES - 1; k++) {
      f[SPECIES_INDEX(k)] = p[SPECIES_INDEX(k)] * p[YVELOCITY];
    }

    for (int i = CONSERVATIVE_VARIABLES; i < CONSERVATIVE_VARIABLES + NONCONSERVATIVE_VARIABLES; i++) {
      f[i]         = p[YVELOCITY] * p[i];
    }
  }
}

template<int blockDimX, int blockDimY, bool X, bool SET>
__global__ void getMUSCLFluxes(Mesh<GPU>::type u, Mesh<GPU>::type flux, const real dt) {
  const bool Y = !X;
  const int dimX = blockDimX - 1 * X,
            dimY = blockDimY - 1 * Y;
  const int i = dimX * blockIdx.x + threadIdx.x - u.ghostCells();
  const int j = dimY * blockIdx.y + threadIdx.y - u.ghostCells();

  const real dx = u.dx();
  const real dy = u.dy();

  __shared__ real temp_shared[NUMBER_VARIABLES][blockDimY][blockDimX];
  real temp1_local[NUMBER_VARIABLES];
  real temp2_local[NUMBER_VARIABLES];
  real temp3_local[NUMBER_VARIABLES];
  Cell temp0(&temp_shared[0][threadIdx.y][threadIdx.x], blockDimX * blockDimY),
       temp1(&temp1_local[0], 1),
       temp2(&temp2_local[0], 1),
       temp3(&temp3_local[0], 1);
  if (u.within(i, j, 0)) {
    Cell& uC = temp0,
        & uL = temp1,
        & uR = temp2;
    // have them all read in the centre point
    for (int k = 0; k < NUMBER_VARIABLES; k++) uC[k] = u(i, j, k);
    __syncthreads();

    if (u.within(i, j, 1)) {

      // read in the left value
      if ((X && threadIdx.x == 0) || (Y && threadIdx.y == 0)) {
        for (int k = 0; k < NUMBER_VARIABLES; k++) uL[k] = u(i - X, j - Y, k);
      } else {
        for (int k = 0; k < NUMBER_VARIABLES; k++) uL[k] = temp_shared[k][threadIdx.y - Y][threadIdx.x - X];
      }

      // read in the right value
      if ((X && threadIdx.x + 1 == blockDimX) || (Y && threadIdx.y + 1 == blockDimY)) {
        for (int k = 0; k < NUMBER_VARIABLES; k++) uR[k] = u(i + X, j + Y, k);
      } else {
        for (int k = 0; k < NUMBER_VARIABLES; k++) uR[k] = temp_shared[k][threadIdx.y + Y][threadIdx.x + X];
      }
      // [1] = u_{i-1}
      // [0] = u_i
      // [2] = u_{i+1}

      // MUSCL reconstruction
      conservativeToPrimitiveInPlace(uL);
      conservativeToPrimitiveInPlace(uC);
      conservativeToPrimitiveInPlace(uR);
      extrapolate(uL, uC, uR);
      primitiveToConservativeInPlace(uL);
      primitiveToConservativeInPlace(uR);
      // [1] = U_L
      // [2] = U_R

      Cell& fL = temp0,
        & fR = temp3;
      // calculate the fluxes---we need to keep U_{L, R} for later, but we no longer need U_i for this calculation.
      conservativeToFlux<X>(uL, fL);
      conservativeToFlux<X>(uR, fR);
      // [0] = F_L
      // [3] = F_R

      // add these fluxes to U_{L,R}
      for (int k = 0; k < NUMBER_VARIABLES; k++) {
        const real delta = dt / (2.0 * (X ? dx : dy))  * (fL[k] - fR[k]);
        uL[k] += delta;
        uR[k] += delta;
      }
      // [1] = \bar{U}_L
      // [2] = \bar{U}_R

      Cell& pL = temp3,
        & pR = temp0;
      // [3] = p_L
      // [0] = p_R

      conservativeToPrimitive(uL, pL);
      conservativeToPrimitive(uR, pR);
      // about to have to use data from adjacent cells, so sync the threads!
      __syncthreads();

      if ((X && threadIdx.x > 0) || (Y && threadIdx.y > 0)) {
       /* if (u(i, j, ISSHOCK)) {
          riemannSolverHLLE<X>(
              Cell(&temp_shared[0][threadIdx.y - Y][threadIdx.x - X], blockDimX * blockDimY),
              temp3,
              0.0, temp1, temp2);
        } else {
*/
          riemannSolverHLLC<X>(
              Cell(&temp_shared[0][threadIdx.y - Y][threadIdx.x - X], blockDimX * blockDimY),
              temp3,
              0.0, temp1, temp2);
       // }

        __syncthreads();

        for (int k = 0; k < NUMBER_VARIABLES; k++) temp0[k] = temp1[k];

        __syncthreads();

        for (int k = 0; k < NUMBER_VARIABLES; k++) {
          if (SET) { flux(i, j, k) = 0.0; }
          flux(i, j, k) += dt / (X ? dx : dy) * temp_shared[k][threadIdx.y][threadIdx.x];
        }
      }
    }
  }
}

template<bool X, bool SET>
__global__ void addFluxesKernel(Mesh<GPU>::type u, Mesh<GPU>::type flux) {
  const bool Y = !X;
  const int i = blockDim.x * blockIdx.x + threadIdx.x;
  const int j = blockDim.y * blockIdx.y + threadIdx.y;

  if (
    (X && i >= -u.ghostCells() + 1 && i < u.activeNx() + u.ghostCells() - 1 && j >= -u.ghostCells() && j < u.activeNy() + u.ghostCells()) ||
    (Y && i >= -u.ghostCells() && i < u.activeNx() + u.ghostCells() && j >= -u.ghostCells() + 1 && j < u.activeNy() + u.ghostCells() - 1)
     ) {	 
#ifdef GHOST
      if (u(i, j, PHI) < 0.0)
#endif
      {
        #pragma unroll
        for (int k = 0; k < CONSERVATIVE_VARIABLES + NONCONSERVATIVE_VARIABLES; k++) {
          if (SET) u(i, j, k) = 0.0;
          u(i, j, k) += flux(i, j, k) - flux(i + X, j + Y, k);
        }
		
	real p[NUMBER_VARIABLES];
    conservativeToPrimitive(u(i, j), p);
      if (p[PRESSURE] > u(i, j, PMAX)) {
         u(i, j, PMAX) = p[PRESSURE];
      }
	  
    }
  }
}

__global__ void addSemiFluxesKernel(Mesh<GPU>::type u, Mesh<GPU>::type flux) {
  const int i = blockDim.x * blockIdx.x + threadIdx.x;
  const int j = blockDim.y * blockIdx.y + threadIdx.y;

  if (u.active(i, j)) {
#ifdef GHOST
    if (u(i, j)[PHI] < 0)
#endif
    {
      #pragma unroll
      for (int k = 0; k < CONSERVATIVE_VARIABLES + NONCONSERVATIVE_VARIABLES; k++) {
        u(i, j, k) += flux(i, j, k);
      }
    }
  }
}

