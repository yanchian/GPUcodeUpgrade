/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "muscl.hpp"

template<int DIRECTION>
__device__ __host__ void starState(const Cell p, const real SS, const real SD, Cell star) {
	const bool X = DIRECTION == XDIRECTION;
	const bool Y = DIRECTION == YDIRECTION;
	const bool Z = DIRECTION == ZDIRECTION;

	int MOMENTUM, TMOMENTUM1, TMOMENTUM2;
	if (X) {
		MOMENTUM   = XMOMENTUM;
		TMOMENTUM1 = YMOMENTUM;
		TMOMENTUM2 = ZMOMENTUM;
	} else if (Y) {
		MOMENTUM   = YMOMENTUM;
		TMOMENTUM1 = XMOMENTUM;
		TMOMENTUM2 = ZMOMENTUM;
	} else if (Z) {
		MOMENTUM   = ZMOMENTUM;
		TMOMENTUM1 = XMOMENTUM;
		TMOMENTUM2 = YMOMENTUM;
	}
	const real chi = (SD - p[MOMENTUM]) / (SD - SS);
	star[DENSITY] = chi * p[DENSITY];
	star[MOMENTUM] = chi * p[DENSITY] * SS;
	star[TMOMENTUM1] = chi * p[DENSITY] * p[TMOMENTUM1];
	star[TMOMENTUM2] = chi * p[DENSITY] * p[TMOMENTUM2];
	star[ENERGY] = chi * ((p[PRESSURE] / (gamma(p) - 1.0) + 0.5 * p[DENSITY] * (p[XVELOCITY] * p[XVELOCITY] + p[YVELOCITY] * p[YVELOCITY] + p[ZVELOCITY] * p[ZVELOCITY])) + (SS - p[MOMENTUM]) * (p[DENSITY] * SS + p[PRESSURE] / (SD - p[MOMENTUM])));
	star[FRACTION] = chi * p[FRACTION];
}

template<int DIRECTION>
__device__ __host__ void riemannSolver(const Cell p_L, const Cell p_R, const real S, Cell solution, Cell temp) {
	const bool X = DIRECTION == XDIRECTION;
	const bool Y = DIRECTION == YDIRECTION;
	//const bool Z = DIRECTION == ZDIRECTION;
	const int VELOCITY = X ? XVELOCITY : (Y ? YVELOCITY : ZVELOCITY);

	const real cL = soundSpeedPrimitive(p_L), cR = soundSpeedPrimitive(p_R);
	const real uROE = (sqrt(p_L[DENSITY]) * p_L[VELOCITY] + sqrt(p_R[DENSITY]) * p_R[VELOCITY]) / (sqrt(p_L[DENSITY]) + sqrt(p_R[DENSITY]));
	const real cROE = (sqrt(p_L[DENSITY]) * cL + sqrt(p_R[DENSITY]) * cR) / (sqrt(p_L[DENSITY]) + sqrt(p_R[DENSITY]));

	const real SL = min(uROE - cROE, p_L[VELOCITY] - cL);
	const real SR = max(uROE + cROE, p_R[VELOCITY] + cR);

	const real SS = (p_R[PRESSURE] - p_L[PRESSURE] + p_L[DENSITY] * p_L[VELOCITY] * (SL - p_L[VELOCITY]) - p_R[DENSITY] * p_R[VELOCITY] * (SR - p_R[VELOCITY])) / (p_L[DENSITY] * (SL - p_L[VELOCITY]) - p_R[DENSITY] * (SR - p_R[VELOCITY]));
	if (SL >= S) {
		primitiveToFlux<DIRECTION>(p_L, solution);
	} else if (SS >= S) {
		primitiveToFlux<DIRECTION>(p_L, solution);
		starState<DIRECTION>(p_L, SS, SL, temp);
		primitiveToConservativeInPlace(p_L);
		for (int k = 0; k < p_L.length; k++) {
			solution[k] += SL * (temp[k] - p_L[k]);
		}
	} else if (SR >= S) {
		primitiveToFlux<DIRECTION>(p_R, solution);
		starState<DIRECTION>(p_R, SS, SR, temp);
		primitiveToConservativeInPlace(p_R);
		for (int k = 0; k < p_R.length; k++) {
			solution[k] += SR * (temp[k] - p_R[k]);
		}
	} else {
		primitiveToFlux<DIRECTION>(p_R, solution);
	}
}

