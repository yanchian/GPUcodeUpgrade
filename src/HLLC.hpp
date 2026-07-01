/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#pragma once

template<bool X>
__device__ __host__ void riemannSolverHLLC(const Cell p_L, const Cell p_R, const real S, Cell solution, Cell temp);

template<bool X>
__device__ __host__ void riemannStateHLLC(const Cell p_L, const Cell p_R, const real S, Cell solution);

