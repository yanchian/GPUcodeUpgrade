/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "boundaryconditions.hpp"

template<BoundaryConditions BCs, bool XDIR, bool downstream>
__global__ void setBoundaryConditionsKernel(Mesh<GPU>::type u) {
  const bool YDIR = !XDIR, upstream = !downstream;

  const int i = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();
  

  if (XDIR && u.exists(i, 0)) {

    const real x = u.x(i);

    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      
      if (downstream) {
    //    u(i, -2, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, 1, k);
    //    u(i, -1, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, 0, k); 
        u(i, -2, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, 1, k);
        u(i, -1, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, 0, k); 
      } else {
        u(i, u.activeNy() + 1, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, u.activeNy() - 2, k);
        u(i, u.activeNy()    , k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(i, u.activeNy() - 1, k);
      }

    }
  
    if  (upstream) {
     // const real U_MAX = M_NUM * sqrt( GAMMA * P_0 / RHO_0 );
      // const real U_LID = - 4.0 * U_MAX * x * x + U_MAX;
      const real U_LID = max_lid_vel;      
      u(i, u.activeNy(), XVELOCITY) = 2.0 * U_LID - u(i, u.activeNy() - 1, XVELOCITY);
      u(i, u.activeNy() + 1, XVELOCITY) = 2.0 * u(i, -u.activeNy(), XVELOCITY) - u(i, u.activeNy() - 1, XVELOCITY);
    }

  }
////////////////////////////////////////////////////////////////

  if (YDIR && u.exists(0, i)) {

    for (int k = 0; k < NUMBER_VARIABLES; k++) {
   
      if (downstream) {
        u(-2, i, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(1, i, k);
        u(-1, i, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(0, i, k); 
      } else {
        u(u.activeNx() + 1, i, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(u.activeNx() - 2, i, k);
        u(u.activeNx()    , i, k) = (k == YMOMENTUM || k == XMOMENTUM ? -1.0 : 1.0) * u(u.activeNx() - 1, i, k);
      }

    }

  }

////////////////////////////////////////////////////////////////
}

/*
template<bool XDIR>
__global__ void setSpecialBoundaryConditionsKernel(Mesh<GPU>::type u) {
  const bool YDIR = !XDIR;

  const int k = blockIdx.x * blockDim.x + threadIdx.x - u.ghostCells();

  if (XDIR && u.exists(k, 0) && k < u.i(XCORNER)) {
    const int j = u.j(YCORNER);
    for (int n = 0; n < 2; n++) {
      u(k, j + n + 1) = u(k, j - n);
      u(k, j + n + 1, YMOMENTUM) = -u(k, j + n + 1, YMOMENTUM);

    }
  }
  else if (!XDIR && u.exists(0, k) && k > u.j(YCORNER)) {
    const int i = u.i(XCORNER);
    for (int n = 0; n < 2; n++) {
      u(i - n - 1, k) = u(i + n, k);
      u(i - n - 1, k, YMOMENTUM) = -u(i - n - 1, k, YMOMENTUM);
    }
  }
}
*/
