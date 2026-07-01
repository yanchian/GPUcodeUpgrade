/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include <math.h>

template<bool DIRECTION>
__host__ __device__ void riemannSolverHLLE(const Cell pL, const Cell pR, const real S, Cell solution, Cell temp);

