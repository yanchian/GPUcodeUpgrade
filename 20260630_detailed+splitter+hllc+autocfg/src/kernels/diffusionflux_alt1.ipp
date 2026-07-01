/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
//#define GHOST

#include <iostream>
#include <iomanip>
#include <stack>
#include <utility>
#include <libconfig.h++>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#ifdef GHOST
#include "ghostfluid.hpp"
#endif
/*
#ifdef REACTIVE
#include "source.hpp"
#include "shockdetect.hpp"
#endif
*/
#include "initialconditions.hpp"
#include "SDF/BoundaryMesh.hpp"
#include "sdf.hpp"

#ifdef GHOST
class Geometry {
private:
  LevelSet<GPU>::type* SDF_;
public:
  bool rotating;
  real omega;

  bool inflow;
  real entropy;
  real velocity;

  LevelSet<GPU>::type& SDF() { return *SDF_; }
  LevelSet<GPU>::type SDF() const { return *SDF_; }
  LevelSet<GPU>::type*& SDFPointer() { return SDF_; }
};

template<typename Mesh>
void SDFShape(Mesh& sdf, const Shape* const shape) {
  for (int i = -sdf.ghostCells(); i < sdf.activeNx() + sdf.ghostCells(); i++) {
    for (int j = -sdf.ghostCells(); j < sdf.activeNy() + sdf.ghostCells(); j++) {
      shape->distance(sdf.x(i, j), sdf(i, j), sdf.d());
    }
  }
}

template<typename Mesh>
void SDFUnary(Mesh& sdf, const SDFUnaryOperation* const op) {
  for (int i = -sdf.ghostCells(); i < sdf.activeNx() + sdf.ghostCells(); i++) {
    for (int j = -sdf.ghostCells(); j < sdf.activeNy() + sdf.ghostCells(); j++) {
      (*op)(sdf(i, j), sdf(i, j));
    }
  }
}
template<typename Mesh>
void SDFBinary(Mesh& sdfa, Mesh& sdfb, const SDFBinaryOperation* const op) {
  for (int i = -sdfa.ghostCells(); i < sdfa.activeNx() + sdfa.ghostCells(); i++) {
    for (int j = -sdfa.ghostCells(); j < sdfa.activeNy() + sdfa.ghostCells(); j++) {
      (*op)(sdfa(i, j), sdfb(i, j), sdfa(i, j));
    }
  }
}
#endif
//////////////////////////////////////////////////////////////////////////////////////
class Solver {
public:
  double timeFluxes, timeDiffusionFluxes, timeSourcing, timeReducing, timeAdding;
  int stepNumber;
  Mesh<GPU>::type* u;
  Mesh<GPU>::type* fluxes;
#ifdef GHOST
  std::vector<Geometry> geometries;
#endif
  std::vector<double> outputTimes;
  std::string outputDirectory;
  std::vector<double> outputRadii;
  real targetCFL;
  int outputNumber;
  double lastDt;
  double omega;

public:
  enum Status {OK, OUTPUT, FINISHED};

private:
  Solver::Status status;

public:
  Solver(std::string filename) :
    lastDt(0.0),
    timeFluxes(0.0),
    timeDiffusionFluxes(0.0),
    timeSourcing(0.0),
    timeReducing(0.0),
    timeAdding(0.0),
    stepNumber(0),
    outputNumber(0)
  {
    using namespace libconfig;

    Config config;
    // automatically cast float <-> int
    config.setAutoConvert(true);
    config.readFile(filename.c_str());

    Setting& root = config.getRoot();

    // Parse simulation parameters
    const Setting& simulation = root["simulation"];

    if (simulation.exists("outputDirectory")) {
      outputDirectory = simulation["outputDirectory"].c_str();
    } else {
      outputDirectory = "";
    }

    int Nx, Ny;
    real xMin, xMax, yMin, yMax;

    Nx = simulation["grid"]["cells"]["x"];
    Ny = simulation["grid"]["cells"]["y"];

    xMin = simulation["grid"]["size"]["x"][0];
    xMax = simulation["grid"]["size"]["x"][1];
    yMin = simulation["grid"]["size"]["y"][0];
    yMax = simulation["grid"]["size"]["y"][1];

    targetCFL = simulation["targetCFL"];

    u = new Mesh<GPU>::type(Nx, Ny, 2, xMin, xMax, yMin, yMax);
    fluxes = new Mesh<GPU>::type(*u, Mesh<GPU>::type::Allocate);

    if (simulation.exists("start")) {
      real start = simulation["start"];
      real end = simulation["end"];
      real interval = simulation["interval"];
      for (int i = 0; start + i * interval < end; i++) {
        outputTimes.push_back(start + i * interval);
      }
    } else {
      outputTimes.push_back(1e30);
    }

    dim3 blockDim(64);
    dim3 gridDim = u->totalCover(blockDim);
    setInitialConditions<<<gridDim, blockDim>>>(*u);

    checkForError();
#ifdef GHOST
    updateLevelSet();
#endif
  }

  Solver::Status step();

  double setBoundaryConditions();

  //template<bool XDIRECTION>
  // double setSpecialBoundaryConditions(void);

  double getDt(real&);
  double getXFluxes(real dt);
  double getYFluxes(real dt);
  double getXDiffusionFluxes(real dt);
  double getYDiffusionFluxes(real dt);
  double addXFluxes();
  double addYFluxes();
  double checkValidity();
#ifdef GHOST
  double updateLevelSet();
  double getTorque(const Geometry& geometry);
  double ghost(bool reflect = true);
  std::pair<real, real> getPressureIntegral(const real radius);
#endif
/*
#ifdef REACTIVE
  double source(real dt);
  void shockDetect();
#endif
*/
  template<bool T>
    double addFluxes();
  double getTimeFluxes() { return timeFluxes; }
  double getTimeDiffusionFluxes() { return timeDiffusionFluxes; }
  double getTimeSourcing() { return timeSourcing; }
  double getTimeReducing() { return timeReducing; }
  double getTimeAdding() { return timeAdding; }
  double getLastDt() { return lastDt; }
  int getStepNumber() { return stepNumber; }
  int getOutputNumber() { return outputNumber; }
};

Solver::Status Solver::step() {
  float timeFlux = 0.0, timeDiffusionFlux = 0.0, timeReduce = 0.0, timeGhost = 0.0, timeBCs = 0.0;

  if (u->time() >= outputTimes[outputTimes.size() - 1]) {
    return FINISHED;
  }

  real dt;
  timeBCs = setBoundaryConditions();

  timeReduce = getDt(dt);

  timeFlux += getXFluxes(dt);
  timeFlux += addXFluxes();
  timeBCs += setBoundaryConditions();

  timeDiffusionFlux += getXDiffusionFluxes(dt);
  timeDiffusionFlux += addXFluxes();
  timeBCs += setBoundaryConditions();

  timeFlux += getYFluxes(dt);
  timeFlux += addYFluxes();
  timeBCs += setBoundaryConditions();

  timeDiffusionFlux += getYDiffusionFluxes(dt);
  timeDiffusionFlux += addYFluxes();
 // timeBCs = setBoundaryConditions();

/*
#ifdef REACTIVE
  source(dt);
#endif
*/
  checkValidity();

  u->time(u->time() + dt);

  timeFluxes += timeFlux;
  timeDiffusionFluxes += timeDiffusionFlux;
  timeAdding += timeBCs;
  timeReducing += timeReduce;
  timeSourcing += timeGhost;

  size_t freeMemory, totalMemory;
  cudaMemGetInfo(&freeMemory, &totalMemory);

  stepNumber++;
  lastDt = dt;
  std::cout << "# Step " << stepNumber << "[" << outputNumber << "]: " << std::fixed << std::setprecision(11) << u->time() << " (dt=" << std::setprecision(3) << std::scientific << dt << "). Time {fluxes=" << std::fixed << std::setprecision(3) << timeFlux << "ms, diffusion fluxes =" << std::fixed << std::setprecision(3) << timeDiffusionFlux << "ms, ghost=" << timeGhost << "ms, reduction=" << timeReduce << "ms, BCs=" << timeBCs << "ms}. Memory " << std::fixed << std::setprecision(0) << (float) ((totalMemory - freeMemory) >> 20) << " MiB / " << (float) (totalMemory >> 20) << " MiB" << std::endl;

  return this->status;
}

double Solver::setBoundaryConditions(void) {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(256);
  dim3 gridDim((max(u->Nx(), u->Ny()) + 255) / 256);
 // dim3 gridDimX((u->Nx() + 255) / 255); 
  setBoundaryConditionsKernel<TRANSMISSIVE, YDIRECTION, DOWNSTREAM><<<gridDim, blockDim>>>(*u);
  cudaThreadSynchronize();
  setBoundaryConditionsKernel<TRANSMISSIVE, YDIRECTION, UPSTREAM><<<gridDim, blockDim>>>(*u);
  cudaThreadSynchronize();

 // dim3 gridDimY((u->Ny() + 255) / 255); 
  setBoundaryConditionsKernel<REFLECTIVE, XDIRECTION, DOWNSTREAM><<<gridDim, blockDim>>>(*u);
  cudaThreadSynchronize();
  setBoundaryConditionsKernel<REFLECTIVE, XDIRECTION, UPSTREAM><<<gridDim, blockDim>>>(*u);
  cudaThreadSynchronize();

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}

/*
template<bool XDIRECTION>
double Solver::setSpecialBoundaryConditions(void) {
  cudaEvent_t start, stop;
  float time;

  checkForError();
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(256);
  dim3 gridDim((max(u->Nx(), u->Ny()) + 255) / 256);
  setSpecialBoundaryConditionsKernel<XDIRECTION><<<gridDim, blockDim>>>(*u);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}
*/

double Solver::getDt(real& dt) {
  cudaEvent_t start, stop;
  float time1;
  float time2;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

	dim3 blockDimReduce(16, 16);
	dim3 gridDimReduce = u->activeCover(blockDimReduce);

        // Find fastest signal speed
	Grid2D<real, GPU, 1> output(gridDimReduce.x, gridDimReduce.y);
	checkForError();

        cudaEventRecord(start, 0);
	getFastestWaveSpeed<256><<<gridDimReduce, blockDimReduce>>>(*u, output);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	checkForError();
	cudaEventElapsedTime(&time1, start, stop);

	Grid2D<real, CPU, 1> outputCPU(output);
	output.free();
	real maximum = 0.0;
	for (int i = 0; i < outputCPU.Nx(); i++) {
	  for (int j = 0; j < outputCPU.Ny(); j++) {
  		maximum = max(maximum, outputCPU(i, j, 0));
          }
	}
	outputCPU.free();

        // Find minimum density
	Grid2D<real, GPU, 1> output1(gridDimReduce.x, gridDimReduce.y);
	checkForError();

        cudaEventRecord(start, 0);
	getMin<256><<<gridDimReduce, blockDimReduce>>>(*u, output1, DENSITY);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	checkForError();
	cudaEventElapsedTime(&time2, start, stop);

	Grid2D<real, CPU, 1> outputCPU1(output1);
	output1.free();
        real min_rho = 9999999.9;
	for (int i = 0; i < outputCPU1.Nx(); i++) {
	  for (int j = 0; j < outputCPU1.Ny(); j++) {
  		min_rho = min(min_rho, outputCPU1(i, j, 0));
          }
	}
	outputCPU1.free();

  real cfl = targetCFL;
  if (stepNumber < 5) {
    cfl *= 0.2;
  }

  real dt_hydro = targetCFL * min(u->dx(), u->dy()) / maximum;
  real dt_diffu = 0.0;

  if ( PR_NUM < 1.0 ) {
    dt_diffu = 0.4 * ( min(u->dx(), u->dy()) * min(u->dx(), u->dy()) ) / ( 2.0 * (mu() / PR_NUM) / min_rho );
  } else {
    dt_diffu = 0.4 * ( min(u->dx(), u->dy()) * min(u->dx(), u->dy()) ) / ( 2.0 * mu() / min_rho );
  }
 // printf("%f %f %f\n", dt_hydro, dt_diffu);
  dt = min(dt_hydro, dt_diffu);
 // dt = dt_hydro;
  if (u->time() + dt >= outputTimes[outputNumber]) {
    dt = outputTimes[outputNumber] - u->time();
    outputNumber++;
    //if (outputNumber == outputTimes.size()) {
    //  status = FINISHED;
    //} else {
    status = OUTPUT;
    //}
  } else {
    status = OK;
  }

  return time1 + time2;
}

// Compute hydrodynamic fluxes /////////////////////////////////////////////////////
double Solver::getXFluxes(real dt) {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(16, 8);
  dim3 gridDim = u->totalCover(blockDim, 1, 0);

  getMUSCLFluxes<16, 8, true, true><<<gridDim, blockDim>>>(*u, *fluxes, dt);
  cudaThreadSynchronize();

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}
double Solver::getYFluxes(real dt) {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(8, 16);
  dim3 gridDim = u->totalCover(blockDim, 0, 1);

  getMUSCLFluxes<8, 16, false, true><<<gridDim, blockDim>>>(*u, *fluxes, dt);
  cudaThreadSynchronize();

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}
////////////////////////////////////////////////////////////////////////////////////

// Compute diffusive fluxes /////////////////////////////////////////////////////
double Solver::getXDiffusionFluxes(real dt) {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(16, 8);
  dim3 gridDim = u->totalCover(blockDim, 1, 0);

  getDiffusionFluxesKernel<16, 8, true, true><<<gridDim, blockDim>>>(*u, *fluxes, dt);
  cudaThreadSynchronize();

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}
double Solver::getYDiffusionFluxes(real dt) {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(8, 16);
  dim3 gridDim = u->totalCover(blockDim, 0, 1);

  getDiffusionFluxesKernel<8, 16, false, true><<<gridDim, blockDim>>>(*u, *fluxes, dt);
  cudaThreadSynchronize();

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}
////////////////////////////////////////////////////////////////////////////////////
double Solver::addXFluxes(void) {
  return addFluxes<true>();
}

double Solver::addYFluxes(void) {
  return addFluxes<false>();
}

template<bool X>
double Solver::addFluxes(void) {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(32, 16);
  dim3 gridDim = u->totalCover(blockDim);
  addFluxesKernel<X, false><<<gridDim, blockDim>>>(*u, *fluxes);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}
double Solver::checkValidity(void) {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(32, 16);
  dim3 gridDim = u->totalCover(blockDim);
  //checkValidityKernel<<<gridDim, blockDim>>>(*u);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}

#ifdef GHOST
double Solver::updateLevelSet() {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(16, 16);
  dim3 gridDim = u->totalCover(blockDim);
  for (int i = 0; i < geometries.size(); i++) {
    if (geometries[i].rotating) {
      combineRotatedLevelSet<<<gridDim, blockDim>>>(*u, geometries[i].SDF(), geometries[i].omega * u->time(), i == 0);
    } else {
      combineLevelSet<<<gridDim, blockDim>>>(*u, geometries[i].SDF(), i == 0);
    }
    cudaThreadSynchronize();
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}

double Solver::getTorque(const Geometry& geometry) {
  dim3 blockDim(16, 16);
  dim3 gridDim = u->totalCover(blockDim);

  Grid2D<real, GPU, 1> torqueField(u->activeNx(), u->activeNy());

  ghostTorque<<<gridDim, blockDim>>>(*u, geometry.SDF(), torqueField, geometry.rotating ?  geometry.omega * u->time() : 0.0);
  cudaThreadSynchronize();

  real torque = thrust::reduce(
      thrust::device_ptr<real>(&(torqueField(0, 0, 0))),
      thrust::device_ptr<real>(&(torqueField(torqueField.Nx() - 1, torqueField.Ny() - 1, 0))),
      0.0);

  torqueField.free();
  return torque;
}
std::pair<real, real> Solver::getPressureIntegral(const real radius) {
  dim3 blockDim(16, 16);
  dim3 gridDim = u->totalCover(blockDim);

  Grid2D<real, GPU, 4> pressureField(u->activeNx(), u->activeNy());

  ghostPressureIntegral<<<gridDim, blockDim>>>(*u, pressureField, radius);
  cudaThreadSynchronize();getYFluxes

  real integrals[4];

  for (int i = 0; i < 4; i++) {
    integrals[i] = thrust::reduce(
        thrust::device_ptr<real>(&(pressureField(0, 0, i))),
        thrust::device_ptr<real>(&(pressureField(pressureField.Nx() - 1, pressureField.Ny() - 1, i))),
        0.0);
  }

  pressureField.free();
  //std::cout << "torque " << std::setprecision(10) << integrals[0] << " " << integrals[1] << " " << integrals[2] << " " << integrals[3] << std::endl;
  return std::pair<real, real>(integrals[0] / integrals[1] + integrals[2] / integrals[3], integrals[3] );
}

double Solver::ghost(const bool reflect) {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(16, 16);
  dim3 gridDim = u->totalCover(blockDim);

  *fluxes = *u;
  for (int i = 0; i < 20; i++) {
    ghostIterate<<<gridDim, blockDim>>>(*u, *fluxes);
    cudaThreadSynchronize();
    swap(*u, *fluxes);10
  }

  if (reflect) {
    ghostReflect<<<gridDim, blockDim>>>(*u, omega);
    cudaThreadSynchronize();

    for (int i = 0; i < geometries.size(); i++) {
      if (geometries[i].rotating) {
        ghostAddRotationalVelocity<<<gridDim, blockDim>>>(*u, geometries[i].SDF(), geometries[i].omega * u->time(), geometries[i].omega);
        cudaThreadSynchronize();
      }
      if (geometries[i].inflow) {
        ghostAddInflow<<<gridDim, blockDim>>>(*u, geometries[i].SDF(), geometries[i].entropy, geometries[i].velocity);
        cudaThreadSynchronize();
      }
    }
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}
#endif

/*
#ifdef REACTIVE
double Solver::source(real dt) {
  cudaEvent_t start, stop;
  float time;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);
  dim3 blockDim(16, 16);
  dim3 gridDim = u->totalCover(blockDim);

  sources<16, 16><<<gridDim, blockDim>>>(*u, dt);
  cudaThreadSynchronize();

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  checkForError();
  cudaEventElapsedTime(&time, start, stop);

  return time;
}

void Solver::shockDetect() {
  dim3 blockDim(16, 16);
  dim3 gridDim = u->activeCover(blockDim);

  shockDetectKernel<<<gridDim, blockDim>>>(*u);
  cudaThreadSynchronize();
  //shockBurnKernel<<<gridDim, blockDim>>>(u);
  //cudaThreadSynchronize();
}
#endif
*/

