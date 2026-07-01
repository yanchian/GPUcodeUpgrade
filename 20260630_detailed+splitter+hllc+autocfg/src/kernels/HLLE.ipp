/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include <math.h>

template<bool DIRECTION>
__host__ __device__ void riemannSolverHLLE(const Cell pL, const Cell pR, const real S, Cell solution, Cell temp) {
  const int VELOCITY = DIRECTION ? XVELOCITY : YVELOCITY;

  real asq,qsq,cfsq,cfl,cfr,bp,bm,ct2=0.0,tmp;
  real al,ar;

  real Gamma = gamma(pL);

  Cell& fL = solution;
  Cell& fR = temp;

  real energyL = pL[PRESSURE] / (Gamma - 1.0) + 0.5 * pL[DENSITY] * (pL[XVELOCITY] * pL[XVELOCITY] + pL[YVELOCITY] * pL[YVELOCITY]);
  real energyR = pR[PRESSURE] / (Gamma - 1.0) + 0.5 * pR[DENSITY] * (pR[XVELOCITY] * pR[XVELOCITY] + pR[YVELOCITY] * pR[YVELOCITY]);

  real rootdL = sqrt(pL[DENSITY]);
  real rootdR = sqrt(pR[DENSITY]);
  real roeNorm = 1.0/(rootdL + rootdR);

  real vRoe[2];
  vRoe[0] = (rootdL * pL[XVELOCITY] + rootdR * pR[XVELOCITY])*roeNorm;
  vRoe[1] = (rootdL * pL[YVELOCITY] + rootdR * pR[YVELOCITY])*roeNorm;
  //vRoe[2] = (rootdL*pL[ZVELOCITY] + rootdR*pR[ZVELOCITY])*roeNorm;


  real hRoe = ((energyL + pL[PRESSURE]) / rootdL + (energyR + pR[PRESSURE]) / rootdR) * roeNorm;

  real aRoe = sqrt(Gamma *  max(hRoe - 0.5 * (vRoe[0] * vRoe[0] + vRoe[1] * vRoe[1]), 1e-8));

  asq = Gamma * pL[PRESSURE] / pL[DENSITY];
  qsq = ct2 + asq;
  tmp = ct2 - asq;
  cfsq = 0.5*(qsq + sqrt((tmp*tmp + 4.0*asq*ct2)));
  cfl = sqrt(cfsq);

  asq = Gamma*pR[PRESSURE]/pR[DENSITY];
  qsq = ct2 + asq;
  tmp = ct2 - asq;
  cfsq = 0.5*(qsq + sqrt((tmp*tmp + 4.0*asq*ct2)));
  cfr = sqrt(cfsq);

/* take max/min of Roe eigenvalues and L/R state wave speeds */
  ar = max(aRoe + vRoe[DIRECTION], (pR[VELOCITY] + cfr));
  al = min(aRoe - vRoe[DIRECTION], (pL[VELOCITY] - cfl));

  bp = max(ar, 0.0);
  bm = min(al, 0.0);

/*--- Step 5. ------------------------------------------------------------------
 * Compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}
 */

  fL[DENSITY]  = pL[VELOCITY] * pL[DENSITY] - bm * pL[DENSITY];
  fR[DENSITY]  = pR[VELOCITY] * pR[DENSITY] - bp * pR[DENSITY];

  fL[XMOMENTUM] = pL[XVELOCITY] * pL[DENSITY] * (pL[VELOCITY] - bm);
  fR[XMOMENTUM] = pR[XVELOCITY] * pR[DENSITY] * (pR[VELOCITY] - bp);

  fL[YMOMENTUM] = pL[YVELOCITY] * pL[DENSITY] * (pL[VELOCITY] - bm);
  fR[YMOMENTUM] = pR[YVELOCITY] * pR[DENSITY] * (pR[VELOCITY] - bp);

  //fL[ZMOMENTUM] = pL[ZVELOCITY] * pL[DENSITY] * (pL[XVELOCITY] - bm);
  //fR[ZMOMENTUM] = pR[ZVELOCITY] * pR[DENSITY] * (pR[XVELOCITY] - bp);

  fL[VELOCITY] += pL[PRESSURE];
  fR[VELOCITY] += pR[PRESSURE];

  fL[ENERGY]  = energyL * (pL[VELOCITY] - bm) + pL[PRESSURE] * pL[VELOCITY];
  fR[ENERGY]  = energyR * (pR[VELOCITY] - bp) + pR[PRESSURE] * pR[VELOCITY];

  for (int k = 0; k < NUM_SPECIES - 1; k++) {
    fL[SPECIES_INDEX(k)] = pL[VELOCITY] * pL[SPECIES_INDEX(k)] - bm * pL[SPECIES_INDEX(k)];
    fR[SPECIES_INDEX(k)] = pR[VELOCITY] * pR[SPECIES_INDEX(k)] - bp * pR[SPECIES_INDEX(k)];
  }

/*--- Step 6. ------------------------------------------------------------------
 * Compute the HLLE flux at interface.
 */

  tmp = 0.5 * (bp + bm) / (bp - bm);
  for (int l = 0; l < NUMBER_VARIABLES; l++) {
    fL[l] = 0.5 * (fL[l] + fR[l]) + (fL[l] - fR[l]) * tmp;
  }
}

