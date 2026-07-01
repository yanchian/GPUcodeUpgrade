#include "diffusionflux.hpp"

template<int blockDimX, int blockDimY, bool X, bool SET>
__global__ void getDiffusionFluxesKernel(Mesh<GPU>::type u, Mesh<GPU>::type flux, const real dt) {
  const bool Y = !X;
  const int dimX = blockDimX - 1 * X,
            dimY = blockDimY - 1 * Y;
  const int i = dimX * blockIdx.x + threadIdx.x - u.ghostCells();
  const int j = dimY * blockIdx.y + threadIdx.y - u.ghostCells();

  const real dx = u.dx();
  const real dy = u.dy();

  //__shared__ real temp_shared[NUMBER_VARIABLES][blockDimY][blockDimX];
  real temp0_local[NUMBER_VARIABLES];
  real temp1_local[NUMBER_VARIABLES];
  real temp2_local[NUMBER_VARIABLES];
  real temp3_local[NUMBER_VARIABLES];
  real temp4_local[NUMBER_VARIABLES];
  real temp5_local[NUMBER_VARIABLES];

  Cell temp0(&temp0_local[0], 1),
       temp1(&temp1_local[0], 1),
       temp2(&temp2_local[0], 1),
       temp3(&temp3_local[0], 1),
       temp4(&temp4_local[0], 1),
       temp5(&temp5_local[0], 1);

  if (u.within(i, j, 0)) {

    Cell& u_I = temp0,
        & u_I_down = temp1,
        & u_I_up = temp2,
        & u_II = temp3,
        & u_II_down = temp4,
        & u_II_up = temp5;
  

    if (u.within(i, j, 1)) {

      // read values from global memory
      for (int k = 0; k < NUMBER_VARIABLES; k++) {

        u_I[k] = u(i, j, k);
        u_I_down[k] = u(i - Y, j - X, k);
        u_I_up[k] = u(i + Y, j + X, k);
        u_II[k] = u(i - X, j - Y, k);
        u_II_down[k] = u(i - X - Y, j - Y - X, k);
        u_II_up[k] = u(i + Y - X, j + X - Y, k);

      }
        
      __syncthreads();

      // Compute diffusive fluxex
      conservativeToPrimitiveInPlace(u_I);
      conservativeToPrimitiveInPlace(u_I_down);
      conservativeToPrimitiveInPlace(u_I_up);
      conservativeToPrimitiveInPlace(u_II);
      conservativeToPrimitiveInPlace(u_II_down);
      conservativeToPrimitiveInPlace(u_II_up);
  
      __syncthreads();

      if ((X && threadIdx.x > 0) || (Y && threadIdx.y > 0)) {

        flux(i, j, DENSITY) = 0.0;

        if (X) {

          flux(i, j, XVELOCITY) = - 1.0 * (2.0/3.0) * mu() * ( 2.0 * ( u_I[XVELOCITY] - u_II[XVELOCITY] ) / dx  - ( ( u_I_up[YVELOCITY] +  u_II_up[YVELOCITY] ) - ( u_I_down[YVELOCITY] +  u_II_down[YVELOCITY] ) ) / (4.0 * dy ) );
  
          flux(i, j, YVELOCITY) = - 1.0 * mu() * ( ( u_I[YVELOCITY] - u_II[YVELOCITY] ) / dx + ( ( u_I_up[XVELOCITY] +  u_II_up[XVELOCITY] ) - ( u_I_down[XVELOCITY] +  u_II_down[XVELOCITY] ) ) / (4.0 * dy ) );
     
        } else {

          flux(i, j, XVELOCITY) = - 1.0 * mu() * ( ( u_I[XVELOCITY] - u_II[XVELOCITY] ) / dy + ( ( u_I_up[YVELOCITY] +  u_II_up[YVELOCITY] ) - ( u_I_down[YVELOCITY] +  u_II_down[YVELOCITY] ) ) / (4.0 * dx ) );

          flux(i, j, YVELOCITY) = - 1.0 * (2.0/3.0) * mu() * ( 2.0 * ( u_I[YVELOCITY] - u_II[YVELOCITY] ) / dy - ( ( u_I_up[XVELOCITY] +  u_II_up[XVELOCITY] ) - ( u_I_down[XVELOCITY] +  u_II_down[XVELOCITY] ) ) / (4.0 * dx ) );
       
        }

        __syncthreads();  
       
        flux(i, j, ENERGY) = flux(i, j, XVELOCITY) * ( u_I[XVELOCITY] + u_II[XVELOCITY] ) / 2.0 + flux(i, j, YVELOCITY) * ( u_I[YVELOCITY] + u_II[YVELOCITY] ) / 2.0 - ( GAMMA_DEFAULT * mu() / PR_NUM ) * ( ( u_I[ENERGY] / (GAMMA_DEFAULT - 1.0) / u_I[DENSITY] ) - ( u_II[ENERGY] / (GAMMA_DEFAULT - 1.0) / u_II[DENSITY] ) ) / (X ? dx : dy);
  
        __syncthreads();

        for (int k = 0; k < NUMBER_VARIABLES; k++) {
        //    if (SET) { flux(i, j, k) = 0.0; }
            flux(i, j, k) *= dt / (X ? dx : dy);
        }

      }
    }
  }
}

