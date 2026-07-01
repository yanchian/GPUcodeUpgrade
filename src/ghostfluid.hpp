#include <cmath>
#include "kernels/SDF.ipp"
#include "Bezier.hpp"
#include "Vector.hpp"
#include "NACA.hpp"

__host__ __device__ real impeller(const Vec p);

__host__ __device__ Vec sdfNormal(const Mesh<GPU>::type& u, const int i, const int j);

__global__ void initialiseDynamicSDF(Mesh<GPU>::type u, Grid2D<real, GPU> dynamicSDF);

__global__ void initialiseStaticSDF(Mesh<GPU>::type u, Grid2D<real, GPU> staticSDF);

__global__ void initialiseLevelSet(Mesh<GPU>::type u, Grid2D<real, GPU> staticSDF, Grid2D<real, GPU> dynamicSDF, real theta);

__host__ __device__ real smoothedSignum(const real x, const real s);

__global__ void ghostTorque(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, Grid2D<real, GPU, 1> integral, const real theta);

__global__ void ghostPressureIntegral(Mesh<GPU>::type u, Grid2D<real, GPU, 4> integral, const real radius);

__global__ void ghostIterate(Mesh<GPU>::type u, Mesh<GPU>::type unext);

__global__ void ghostReflect(Mesh<GPU>::type u, real omega);

__global__ void ghostAddRotationalVelocity(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, real theta, real omega);

__global__ void ghostAddTranslationalVelocity(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, Vec offset, Vec velocity);

__global__ void ghostAddInflow(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, real entropy, real velocity);

__global__ void combineLevelSet(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, bool first);

__global__ void combineRotatedLevelSet(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, real theta, bool first);

__global__ void combineTranslatedLevelSet(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, Vec offset, bool first);

template<typename Mesh>
void fillFractions(Mesh& sdf, Polygon<real>& polygon);

__global__ void checkValidityKernel(Mesh<GPU>::type u);

