/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#include "ghostfluid.hpp"


__host__ __device__ Vec sdfNormal(const Mesh<GPU>::type& u, const int i, const int j) {
  Vec n;
  n[0] = (u(i + 1, j)[PHI] - u(i - 1, j)[PHI]) / (2.0 * u.dx());
  n[1] = (u(i, j + 1)[PHI] - u(i, j - 1)[PHI]) / (2.0 * u.dy());
  return n / abs(n);
}

__host__ __device__ real smoothedDelta(const real x, const real s) {
  return 1.0 / sqrt(2.0 * M_PI) * (1.0 / s) * exp(-0.5 * pow(x / s, 2.0));
}

__host__ __device__ real smoothedSignum(const real x, const real s) {
  return x / sqrt(x * x + s * s);
}

__host__ __device__ void bilinearInterpolate(const Mesh<GPU>::type& u, const Vec x, Cell v) {
  real ri = u.i(x[0]), rj = u.j(x[1]);
  real oi = ri - (int) ri, oj = rj - (int) rj;

  v = 0.0;
  if (u.active(ri, rj) && u.active(ri + 1, rj + 1)) {
    for (int m = 0; m < NUMBER_VARIABLES; m++) {
      v[m] = u(ri, rj, m) * (1.0 - oi) * (1.0 - oj) +
             u(ri + 1, rj, m) * oi * (1.0 - oj) +
             u(ri, rj + 1, m) * (1.0 - oi) * oj +
             u(ri + 1, rj + 1, m) * oi * oj;
    }
  }
}

__global__ void ghostTorque(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, Grid2D<real, GPU, 1> integral, const real theta) {
  int i = blockDim.x * blockIdx.x + threadIdx.x - u.ghostCells();
  int j = blockDim.y * blockIdx.y + threadIdx.y - u.ghostCells();

  if (u.within(i, j, 3)) {
    Vec x = u.x(i, j);

    Vec n = sdfNormal(u, i, j);

    Vec r = rotateCoordinates(x, Vec(0.0, 0.0), -theta);
    real ri = sdf.i(r[0]), rj = sdf.j(r[1]);
    real oi = ri - (int) ri, oj = rj - (int) rj;
    real value;
    if (sdf.active(ri, rj) && sdf.active(ri + 1, rj + 1)) {
      value = sdf(ri, rj, 0) * (1.0 - oi) * (1.0 - oj) +
              sdf(ri + 1, rj, 0) * oi * (1.0 - oj) +
              sdf(ri, rj + 1, 0) * (1.0 - oi) * oj +
              sdf(ri + 1, rj + 1, 0) * oi * oj;
    } else {
      value = -1e30;
    }

    real temp[NUMBER_VARIABLES];
    Cell p(temp);

    conservativeToPrimitive(u(i, j), p);
    integral(i, j, 0) = p[PRESSURE] * (x[0] * n[1] - x[1] * n[0]) * smoothedDelta(value/* + 0.5 * u.dx()*/, u.dx()) * u.dx() * u.dy();
    if (isnan(integral(i, j, 0))) {
      //printf("%d, %d: pressure: %g, x: %g, %g, n: %g, %g, v: %g, int: %g\n", i, j, p[PRESSURE], x[0], x[1], n[0], n[1], value, integral(i, j, 0));
    }
  } else if (u.exists(i, j)) {
    integral(i, j, 0) = 0.0;
  }
}

__global__ void ghostPressureIntegral(Mesh<GPU>::type u, Grid2D<real, GPU, 4> integral, const real radius) {
  int i = blockDim.x * blockIdx.x + threadIdx.x - u.ghostCells();
  int j = blockDim.y * blockIdx.y + threadIdx.y - u.ghostCells();

  if (u.active(i, j)) {
    Vec x = u.x(i, j);

    Vec n = x / abs(x);

    real temp[NUMBER_VARIABLES];
    Cell p(temp);

    conservativeToPrimitive(u(i, j), p);

    const real dl = (1.0 - u(i, j, FRACTION)) * smoothedDelta(abs(x) - radius, u.dx()) * u.dx() * u.dy();

    if (p[PHI] < 0.0) {
      integral(i, j, 0) = p[PRESSURE] * (1.0 - u(i, j, FRACTION)) * u.dx() * u.dy(); // FIX
      integral(i, j, 1) = dl;
      integral(i, j, 2) = 0.5 * p[DENSITY] * pow(abs(p.velocity()),2) * dot(p.velocity(), n) *  dl;
      integral(i, j, 3) = dot(p.velocity(), n) * dl;
    } else {
      integral(i, j, 0) = 0.0;
      integral(i, j, 1) = 0.0;
      integral(i, j, 2) = 0.0;
      integral(i, j, 3) = 0.0;
    }
  }
}

__global__ void ghostIterate(Mesh<GPU>::type u, Mesh<GPU>::type unext) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  const real epsilon = 3 * u.dx();

  if (u.active(i, j)) {
    if (u(i, j)[PHI] > 0 && u(i, j)[PHI] < epsilon) {
      Vec n = sdfNormal(u, i, j);
      for (int l = 0; l < CONSERVATIVE_VARIABLES + NONCONSERVATIVE_VARIABLES; l++) {
        Vec I;
        if (n[0] > 0) {
          I[0] = (u(i, j, l) - u(i - 1, j, l)) / u.dx();
        } else {
          I[0] = (u(i + 1, j, l) - u(i, j, l)) / u.dx();
        }
        if (n[1] > 0) {
          I[1] = (u(i, j, l) - u(i, j - 1, l)) / u.dy();
        } else {
          I[1] = (u(i, j + 1, l) - u(i, j, l)) / u.dy();
        }
        unext(i, j, l) = u(i, j, l) - u.dx() / 2.0 * dot(n, I);
      }
    }
  }
}

__global__ void ghostReflect(Mesh<GPU>::type u, real omega) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  const real epsilon = 3 * u.dx();

  if (u.active(i, j)) {
    if (u(i, j)[PHI] > 0 && u(i, j)[PHI] < epsilon) {
      Vec x = u.x(i, j);
      Vec velocity(0.0);
      Vec n = sdfNormal(u, i, j);
      Vec momentum = u(i, j).momentum();
      momentum = momentum - 2.0 * u(i, j, FRACTION) * dot(n, momentum) * n;
      u(i, j).momentum() = momentum;
    }
  }
}

__global__ void ghostAddRotationalVelocity(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, real theta, real omega) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  const real epsilon = 3 * u.dx();

  if (u.active(i, j)) {
    Vec x = u.x(i, j);
    Vec r = rotateCoordinates(x, Vec(0.0, 0.0), -theta);
    real ri = sdf.i(r[0]), rj = sdf.j(r[1]);
    real oi = ri - (int) ri, oj = rj - (int) rj;

    real value;
    if (sdf.active(ri, rj) && sdf.active(ri + 1, rj + 1)) {
      value = sdf(ri, rj, 0) * (1.0 - oi) * (1.0 - oj) +
              sdf(ri + 1, rj, 0) * oi * (1.0 - oj) +
              sdf(ri, rj + 1, 0) * (1.0 - oi) * oj +
              sdf(ri + 1, rj + 1, 0) * oi * oj;
    } else {
      value = -1e5;
    }
    if (value > 0 && value < epsilon) {
      Vec x = u.x(i, j);
      Vec velocity(0.0);
      velocity[0] = -x[1] * omega;
      velocity[1] =  x[0] * omega;

      Vec n = sdfNormal(u, i, j);
      Vec momentum = u(i, j).momentum();
      momentum = momentum + 2.0 * u(i, j, FRACTION) * dot(n, velocity) * u(i, j, DENSITY) * n;
      u(i, j).momentum() = momentum;
    }
  }
}

__global__ void ghostAddTranslationalVelocity(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, Vec offset, Vec velocity) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  const real epsilon = 3 * u.dx();

  if (u.active(i, j)) {
    Vec x = u.x(i, j);
    Vec r = x - offset;
    real ri = sdf.i(r[0]), rj = sdf.j(r[1]);
    real oi = ri - (int) ri, oj = rj - (int) rj;

    real value;
    if (sdf.active(ri, rj) && sdf.active(ri + 1, rj + 1)) {
      value = sdf(ri, rj, 0) * (1.0 - oi) * (1.0 - oj) +
              sdf(ri + 1, rj, 0) * oi * (1.0 - oj) +
              sdf(ri, rj + 1, 0) * (1.0 - oi) * oj +
              sdf(ri + 1, rj + 1, 0) * oi * oj;
    } else {
      value = -1e5;
    }
    if (value > 0 && value < epsilon) {
      Vec x = u.x(i, j);

      Vec n = sdfNormal(u, i, j);
      Vec momentum = u(i, j).momentum();
      momentum = momentum + 2.0 * u(i, j, FRACTION) * u(i, j, DENSITY) * dot(n, velocity) * n;
      u(i, j).momentum() = momentum;
    }
  }
}

__global__ void ghostAddInflow(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, real entropy, real velocity) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  const real epsilon = 5 * u.dx();

  if (u.active(i, j)) {
    Vec x = u.x(i, j);
    if (sdf(i, j, 0) > 0 && sdf(i, j, 0) < epsilon) {
      conservativeToPrimitiveInPlace(u(i, j));
      Vec n = sdfNormal(u, i, j);
      Vec t(-n[1], n[0]);
      u(i, j).momentum() = n * velocity + t * dot(u(i, j).momentum(), t);
      u(i, j, PRESSURE) = entropy * pow(u(i, j, DENSITY), gamma(u(i, j))) - p0(u(i, j));
      primitiveToConservativeInPlace(u(i, j));
    }
  }
}

__global__ void combineLevelSet(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, bool first) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  if (u.active(i, j)) {
    if (first || u(i, j, PHI) < sdf(i, j, 0)) {
      u(i, j, PHI) = sdf(i, j, 0);
      u(i, j, FRACTION) = sdf(i, j, 1);
    }
      if (isnan(u(i, j, PHI))) {
        //printf("nan at %d, %d\n", i, j);
      }
  }
}

__global__ void combineRotatedLevelSet(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, real theta, bool first) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  if (u.active(i, j)) {
    if (first) {
      u(i, j, PHI) = -1e30;
    }
    Vec x = u.x(i, j);
    Vec r = rotateCoordinates(x, Vec(0.0, 0.0), -theta);
    real ri = sdf.i(r[0]), rj = sdf.j(r[1]);
    real oi = ri - (int) ri, oj = rj - (int) rj;

    real value, frac;
    if (sdf.active(ri, rj) && sdf.active(ri + 1, rj + 1)) {
      value = sdf(ri, rj, 0) * (1.0 - oi) * (1.0 - oj) +
              sdf(ri + 1, rj, 0) * oi * (1.0 - oj) +
              sdf(ri, rj + 1, 0) * (1.0 - oi) * oj +
              sdf(ri + 1, rj + 1, 0) * oi * oj;
      frac = sdf(ri, rj, 1) * (1.0 - oi) * (1.0 - oj) +
             sdf(ri + 1, rj, 1) * oi * (1.0 - oj) +
             sdf(ri, rj + 1, 1) * (1.0 - oi) * oj +
             sdf(ri + 1, rj + 1, 1) * oi * oj;
    } else {
      value = -1e30;
      frac = 0.0;
    }
    if (first || u(i, j, PHI) < value) {
      u(i, j, PHI) = value;
      u(i, j, FRACTION) = frac;
    }
  }
}

__global__ void combineTranslatedLevelSet(Mesh<GPU>::type u, LevelSet<GPU>::type sdf, Vec offset, bool first) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  if (u.active(i, j)) {
    if (first) {
      u(i, j, PHI) = -1e30;
    }
    Vec x = u.x(i, j);
    Vec r = x - offset;
    real ri = sdf.i(r[0]), rj = sdf.j(r[1]);
    real oi = ri - (int) ri, oj = rj - (int) rj;

    real value, frac;
    if (sdf.active(ri, rj) && sdf.active(ri + 1, rj + 1)) {
      value = sdf(ri, rj, 0) * (1.0 - oi) * (1.0 - oj) +
              sdf(ri + 1, rj, 0) * oi * (1.0 - oj) +
              sdf(ri, rj + 1, 0) * (1.0 - oi) * oj +
              sdf(ri + 1, rj + 1, 0) * oi * oj;
      frac = sdf(ri, rj, 1) * (1.0 - oi) * (1.0 - oj) +
             sdf(ri + 1, rj, 1) * oi * (1.0 - oj) +
             sdf(ri, rj + 1, 1) * (1.0 - oi) * oj +
             sdf(ri + 1, rj + 1, 1) * oi * oj;
    } else {
      value = -1e30;
      frac = 0.0;
    }
    if (first || u(i, j, PHI) < value) {
      u(i, j, PHI) = value;
      u(i, j, FRACTION) = frac;
    }
  }
}

template<typename Mesh>
void fillFractions(Mesh& sdf, Polygon<real>& polygon) {
  for (int i = 0; i < sdf.activeNx(); i++) {
    for (int j = 0; j < sdf.activeNy(); j++) {
      if (std::abs(sdf(i, j, 0)) < 2.0 * sdf.dx()) {
        std::vector<Vec> points;
        Vec lastNormal;

        for (int v = 0; v < polygon.points(); v++) {
          Polygon<real>::Point start = polygon.point(v);
          Polygon<real>::Point end   = polygon.point(v + 1);

          const Vec x = sdf.x(i, j);
          Vec corners[2];
          corners[0] = x - sdf.d() / 2.0;
          corners[1] = x + sdf.d() / 2.0;
          if (corners[0][0] < start[0] && start[0] < corners[1][0] && corners[0][1] < start[1] && start[1] < corners[1][1]) {
            points.push_back(start);
          }
          for (int c1 = 0; c1 < 2; c1++) {
            const int c2 = (c1 + 1) % 2;
            for (int e1 = 0; e1 < 2; e1++) {
              const int e2 = (e1 + 1) % 2;
              if ((start[e1] <= corners[c1][e1]) != (end[e1] <= corners[c1][e1])) {
                const real intersection = (start[e2] * (end[e1] - corners[c1][e1]) + end[e2] * (corners[c1][e1] - start[e1])) / (end[e1] - start[e1]);
                if ((corners[c1][e2] <= intersection) != (corners[c2][e2] <= intersection)) {
                  Vec p;
                  p[e1] = corners[c1][e1];
                  p[e2] = intersection;
                  if (v == polygon.points() - 1) {
                    points.insert(points.begin(), p);
                  } else {
                    points.push_back(p);
                  }
                  lastNormal = end - start;
                  lastNormal = Vec(-lastNormal[1], lastNormal[0]);
                  lastNormal /= abs(lastNormal);
                }
              }
            }
          }
        }
        if (points.size() > 0) {
          Vec first = points[0];
          Vec last = points[points.size() - 1];
          bool inside = false, reverse = false;
          for (int k = 0; k < 8; k++) {
            const int e = k % 4;
            Vec c;
            switch (e) {
              case 0: c = Vec(sdf.x(i) - sdf.dx() / 2.0, sdf.y(j) - sdf.dy() / 2.0); break;
              case 1: c = Vec(sdf.x(i) - sdf.dx() / 2.0, sdf.y(j) + sdf.dy() / 2.0); break;
              case 2: c = Vec(sdf.x(i) + sdf.dx() / 2.0, sdf.y(j) + sdf.dy() / 2.0); break;
              case 3: c = Vec(sdf.x(i) + sdf.dx() / 2.0, sdf.y(j) - sdf.dy() / 2.0); break;
            }
            if (inside) {
              points.push_back(c);
              if (c[0] == first[0] || c[1] == first[1]) {
                break;
              }
            }
            if (!inside && ((c[0] == last[0] && e % 2 == 1) || (c[1] == last[1] && e % 2 == 0))) {
              if (dot(c - last, lastNormal) < 0) {
                reverse = true;
              }
              points.push_back(c);
              inside = true;
            }
          }
          real area = 0.0;
          //std::cout << i << ", " << j << " ";
          for (int v = 0; v < points.size(); v++) {
            //std::cout << std::setprecision(10) << "{" << points[v][0] << ", " << points[v][1] << "}, ";
            real a = points[v][0] * points[(v + 1) % points.size()][1] - points[(v + 1) % points.size()][0] * points[v][1];
            area += a;
          }
          //std::cout << " area " << (std::abs(area) / 2.0) / (sdf.dx() * sdf.dy()) << std::endl;
          sdf(i, j, 1) = (std::abs(area) / 2.0) / (sdf.dx() * sdf.dy());
          if (reverse) {
            sdf(i, j, 1) = 1.0 - sdf(i, j, 1);
          }
        } else {
          sdf(i, j, 1) = sdf(i, j, 0) > 0;
        }
      } else {
        sdf(i, j, 1) = sdf(i, j, 0) > 0;
      }
    }
  }
}

__global__ void checkValidityKernel(Mesh<GPU>::type u) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  int j = blockDim.y * blockIdx.y + threadIdx.y;

  if (u.active(i, j)) {
    for (int k = 0; k < NUMBER_VARIABLES; k++) {
      if (isnan(u(i, j, k))) {
        //printf("Invalid data at (%d, %d, %d)!\n", i, j, k);
      }
    }
  }
}

