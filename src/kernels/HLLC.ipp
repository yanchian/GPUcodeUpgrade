/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "HLLC.hpp"

template<bool X>
__device__ __host__ void starState(const Cell p, const real SS, const real SD, Cell star,
                                   const real v_t_override = -1e30) {
  const int MOMENTUM = X ? XMOMENTUMFLUX : YMOMENTUMFLUX;
  const int TMOMENTUM = X ? YMOMENTUMFLUX : XMOMENTUMFLUX;
  const int VELOCITY = X ? XVELOCITY : YVELOCITY;
  const int TVELOCITY = X ? YVELOCITY : XVELOCITY;
  const real chi = (SD - p[VELOCITY]) / (SD - SS);
  star[DENSITY] = chi * p[DENSITY];
  star[MOMENTUM] = chi * p[DENSITY] * SS;

  // 剪切修正: 使用 Roe 平均切向速度替代单侧切向速度
  // 粘性流中接触间断处速度应连续, 否则粘性应力计算不正确
  real v_t = (v_t_override > -1e29) ? v_t_override : p[TVELOCITY];
  star[TMOMENTUM] = chi * p[DENSITY] * v_t;

  real ke = p[XVELOCITY] * p[XVELOCITY] + p[YVELOCITY] * p[YVELOCITY];
  // 若使用了剪切修正, 需要调整动能中的切向分量
  if (v_t_override > -1e29) {
    real v_t_old = p[TVELOCITY];
    ke += (v_t * v_t - v_t_old * v_t_old);
  }
  star[ENERGY] = chi * (((p[PRESSURE] + gamma(p) * p0(p)) / (gamma(p) - 1.0) + 0.5 * p[DENSITY] * ke) + (SS - p[VELOCITY]) * (p[DENSITY] * SS + p[PRESSURE] / (SD - p[VELOCITY])));
  for (int k = 0; k < NUM_SPECIES - 1; k++) {
    star[SPECIES_INDEX(k)] = chi * p[SPECIES_INDEX(k)];
  }
  for (int i = CONSERVATIVE_VARIABLES; i < CONSERVATIVE_VARIABLES + NONCONSERVATIVE_VARIABLES; i++) {
    star[i] = p[i];
  }
}

template<bool X>
__device__ __host__ void riemannSolverHLLC(const Cell p_L, const Cell p_R, const real S, Cell solution, Cell temp) {
  const int VELOCITY = X ? XVELOCITY : YVELOCITY;
  const int TVELOCITY = X ? YVELOCITY : XVELOCITY;

  const real cL = soundSpeedPrimitive(p_L), cR = soundSpeedPrimitive(p_R);
  const real rootdL = sqrt(p_L[DENSITY]), rootdR = sqrt(p_R[DENSITY]);
  const real uROE = (rootdL * p_L[VELOCITY] + rootdR * p_R[VELOCITY]) / (rootdL + rootdR);
  const real cROE = (rootdL * cL + rootdR * cR) / (rootdL + rootdR);

  // Roe 平均切向速度 (粘性流中接触间断处速度应连续)
  const real v_t_roe = (rootdL * p_L[TVELOCITY] + rootdR * p_R[TVELOCITY]) / (rootdL + rootdR);

  const real SL = min(uROE - cROE, p_L[VELOCITY] - cL);
  const real SR = max(uROE + cROE, p_R[VELOCITY] + cR);

  const real SS = (p_R[PRESSURE] - p_L[PRESSURE] + p_L[DENSITY] * p_L[VELOCITY] * (SL - p_L[VELOCITY]) - p_R[DENSITY] * p_R[VELOCITY] * (SR - p_R[VELOCITY])) / (p_L[DENSITY] * (SL - p_L[VELOCITY]) - p_R[DENSITY] * (SR - p_R[VELOCITY]));
  if (SL >= S) {
    primitiveToFlux<X>(p_L, solution);
  } else if (SS >= S) {
    primitiveToFlux<X>(p_L, solution);
    starState<X>(p_L, SS, SL, temp, v_t_roe);
    for (int k = 0; k < p_L.length(); k++) {
      solution[k] += SL * temp[k];
    }
    primitiveToConservative(p_L, temp);
    for (int k = 0; k < p_L.length(); k++) {
      solution[k] -= SL * temp[k];
    }
  } else if (SR >= S) {
    primitiveToFlux<X>(p_R, solution);
    starState<X>(p_R, SS, SR, temp, v_t_roe);
    for (int k = 0; k < p_R.length(); k++) {
      solution[k] += SR * temp[k];
    }
    primitiveToConservative(p_R, temp);
    for (int k = 0; k < p_R.length(); k++) {
      solution[k] -= SR * temp[k];
    }
  } else {
    primitiveToFlux<X>(p_R, solution);
  }
}

template<bool X>
__device__ __host__ void riemannStateHLLC(const Cell p_L, const Cell p_R, const real S, Cell solution) {
  const int VELOCITY = X ? XVELOCITY : YVELOCITY;
  const int TVELOCITY = X ? YVELOCITY : XVELOCITY;

  const real cL = soundSpeedPrimitive(p_L), cR = soundSpeedPrimitive(p_R);
  const real rootdL = sqrt(p_L[DENSITY]), rootdR = sqrt(p_R[DENSITY]);
  const real uROE = (rootdL * p_L[VELOCITY] + rootdR * p_R[VELOCITY]) / (rootdL + rootdR);
  const real cROE = (rootdL * cL + rootdR * cR) / (rootdL + rootdR);

  const real v_t_roe = (rootdL * p_L[TVELOCITY] + rootdR * p_R[TVELOCITY]) / (rootdL + rootdR);

  const real SL = min(uROE - cROE, p_L[VELOCITY] - cL);
  const real SR = max(uROE + cROE, p_R[VELOCITY] + cR);

  const real SS = (p_R[PRESSURE] - p_L[PRESSURE] + p_L[DENSITY] * p_L[VELOCITY] * (SL - p_L[VELOCITY]) - p_R[DENSITY] * p_R[VELOCITY] * (SR - p_R[VELOCITY])) / (p_L[DENSITY] * (SL - p_L[VELOCITY]) - p_R[DENSITY] * (SR - p_R[VELOCITY]));
  if (SL >= S) {
    solution = p_L;
  } else if (SS >= S) {
    starState<X>(p_L, SS, SL, solution, v_t_roe);
    conservativeToPrimitiveInPlace(solution);
  } else if (SR >= S) {
    starState<X>(p_R, SS, SR, solution, v_t_roe);
    conservativeToPrimitiveInPlace(solution);
  } else {
    solution = p_R;
  }
}

