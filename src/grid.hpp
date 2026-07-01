#pragma once
#include <algorithm>
#include <iostream>
#include <ctime>
#include <cstdio>
#include <string.h>
#include <string>
#include "Vector.hpp"

double typedef real;

Vector<real, 2> typedef Vec;

#ifndef NUM_SPECIES
#define NUM_SPECIES 42
#endif

#ifndef NUMBER_VARIABLES
const int NUMBER_VARIABLES = 4 + (NUM_SPECIES - 1) + 2;  // ρ,ρu,ρv,ρE + (N-1)species + ISSHOCK + PMAX
#endif
const int GHOST_CELLS = 2;
const int CONSERVATIVE_VARIABLES = 4 + (NUM_SPECIES - 1);  // all hydrodynamic + species are conservative
const int NONCONSERVATIVE_VARIABLES = 0;
const int DUMMY_VARIABLES = 0;

// Macro to compute species variable index from species k
#define SPECIES_INDEX(k) (4 + (k))

enum ConservedVariables {
  DENSITY, XMOMENTUM, YMOMENTUM, ENERGY,
  ISSHOCK = NUMBER_VARIABLES - 2,
  PMAX = NUMBER_VARIABLES - 1
  // Species variables: SPECIES_INDEX(0) through SPECIES_INDEX(NUM_SPECIES-2)
  // The last species (NUM_SPECIES-1) is not stored; Y_last = 1 - sum(Y_k)
};
enum FluxVariables {DENSITYFLUX, XMOMENTUMFLUX, YMOMENTUMFLUX, ENERGYFLUX};
enum PrimitiveVariables {DENSITYPRIM, XVELOCITY, YVELOCITY, PRESSURE};
enum BoundaryConditions {REFLECTIVE, TRANSMISSIVE};

#include "Timer.hpp"
#include "StridedArray.hpp"
#include "Grid2.hpp"
#include "Mesh2.hpp"

