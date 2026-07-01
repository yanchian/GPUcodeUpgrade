#pragma once
#include "Vector.hpp"

/*
  CNS Hydrocode Package
  Copyright © Cambridge Numerical Solutions Ltd 2009

  SVN-Id: $Id: Bezier.C 1330 2012-12-04 12:24:44Z philipblakely $
 */

__host__ __device__ real distToLine(Vec a, Vec b, Vec p) {
	p = p - a;
	b = b - a;
	const real pb = dot(p, b);
	const real bb = norm2(b);
	real dist = 0;
	if (pb > bb) {
		p = p - b;
		dist = norm(p);
	} else if (pb < 0) {
		dist = norm(p);
	} else {
		dist = abs(p[0] * b[1] - p[1] * b[0]) / sqrt(bb);
	}
	return ((p[0] * b[1] - p[1] * b[0]) >= 0 ? -1 : 1) * dist;
}

template<int deg>
class Bezier {
  //! Bezier curve's control points
  Vec vertices[deg + 1];

  unsigned int m_binCoeff[deg + 1];
	Vec m_binCoeffV[deg + 1];
  Vec m_binCoeffd1V[deg + 1];
  Vec m_binCoeffd2V[deg + 1];
  Vec m_binCoeffdd1V[deg + 1];
  Vec m_binCoeffdd2V[deg + 1];
  Vec m_binCoeffdd3V[deg + 1];

public:
	__host__ __device__ Bezier(const real params[deg + 1][2]) {
		for (size_t i = 0; i <= deg; i++) {
			vertices[i][0] = params[i][0];
			vertices[i][1] = params[i][1];
		}

		m_binCoeff[0] = 1;
		for (unsigned int i = 0; i <= deg; i++) {
			// Compute (n r) and avoid integer overflow
			if (i > 0) {
				m_binCoeff[i] = ((deg - (i - 1)) * m_binCoeff[i - 1]) / i;
			}
			unsigned int c = m_binCoeff[i];
			m_binCoeffV[i] = c * vertices[i];
			m_binCoeffd1V[i] = c * i * vertices[i];
			m_binCoeffd2V[i] = c * (deg - i) * vertices[i];
			m_binCoeffdd1V[i] = c * i * (i - 1) * vertices[i];
			m_binCoeffdd2V[i] = c * 2 * i * (deg - i) * vertices[i];
			m_binCoeffdd3V[i] = c * (deg - i) * (deg - i - 1) * vertices[i];
		}
	}

	__host__ __device__ real distance(const Vec p) const {
		// First we get a good starting point for our iterative method by
		// approximating the Bezier curve as a series of straight lines and finding the closest
		// point to (x,y) on it.

		real minDist = 1e15;
		real nearestW = -1;
		unsigned int parts = 30;
		real h = 1.0 / parts;

		for (unsigned int i = 0 ; i <= parts ; i++) {
			real u = h * i;
			Vec curveu;
			curveNoDerivs(u, curveu);
			const real d = norm(p - curveu);
			if (abs(d) < minDist) {
				minDist = abs(d);
				nearestW = u;
			}
		}

		// Now nearestW is the proportion of the way round the polygon
		// of the nearest point on the polygon to (x,y)

		real u = min(max(nearestW, 0.), 1.);
		real oldU;
		Vec c, dc, d2c, cmp;
		bool endCondn;
		unsigned int iter = 0;
		do {
			curve(u, c, dc, d2c);
			cmp = c - p;
			real f = dot(dc, cmp);
			real df = dot(d2c, cmp) + norm2(dc);

			oldU = u;
			u = oldU - f / df;

			if (u < 1e-6 || u > 1 - 1e-6) {
				break;
			}

			endCondn = (u - oldU) * norm(dc) < 1e-4;
		 	endCondn = endCondn || norm(c - p) < 1e-4;
			endCondn = endCondn || abs(f) / (norm(dc) * norm(cmp)) < 1e-4;
			iter++;
		} while (!endCondn && iter < 10);

		// If we go off the end of the curve,
		// then use the tangent at the end of the curve
		// to determine distance from it, and whether we're above or below.

		const bool offAtLeft = (u < 1e-6);
		const bool offAtRight = (u > 1 -1e-6);
		real offEndDist = 1e15; // inf
		if (offAtLeft || offAtRight) {
			real d1 = distToLine(vertices[0], vertices[1], p);
			real d2 = distToLine(vertices[deg - 1], vertices[deg], p);

			// We use the left-side if we're closest to the beginning, otherwise use right-side
			// Can't depend on value of u, iteration procedure has gone wrong at this stage and
			// could have gone off at left when in fact right-most point is closer.
			if (norm(vertices[0] - p) < norm(vertices[deg] - p)) {
				offEndDist = d1;
			} else {
				offEndDist = d2;
			}
		}

		// Otherwise use the nearest point on curve as found above
		// and use the normal to the curve to determine whether we're above or below.
		u = min(max(u, 0.), 1.);

		curve(u, c, dc, d2c);
		real curveDist;
		if (dc[1] * (c[0] - p[0]) - dc[0] * (c[1] - p[1]) > 0) {
			curveDist = norm(c - p);
		} else {
			curveDist = -norm(c - p);
		}

		if (abs(curveDist) < abs(offEndDist)) {
			return curveDist;
		} else {
			return offEndDist;
		}
	}

//! Evaluate the Bezier curve for parameter u
	__host__ __device__ void curveNoDerivs(real u, Vec& pos) const {
		real uPow[deg + 1];
		real uM1Pow[deg + 1];
		uPow[0] = 1.0;
		uM1Pow[0] = 1.0;
		for (unsigned int i = 1; i <= deg; i++) {
			uPow[i] = uPow[i - 1] * u;
			uM1Pow[i] = uM1Pow[i - 1] * (1 - u);
		}

		pos[0] = pos[1] = 0.0;
		for (size_t i = 0; i <= deg; i++) {
			pos += m_binCoeffV[i] * (uPow[i] * uM1Pow[deg - i]);
		}
	}

//! Evaluate the Bezier curve for parameter u as well as its first and second derivatives
	__host__ __device__ void curve(real u, Vec& pos, Vec& dpos, Vec& d2pos) const {
		real uPow[deg + 1];
		real uM1Pow[deg + 1];
		uPow[0] = 1.0;
		uM1Pow[0] = 1.0;
		for (unsigned int i = 1; i <= deg; i++) {
			uPow[i] = uPow[i - 1] * u;
			uM1Pow[i] = uM1Pow[i - 1] * (1 - u);
		}

		pos = dpos = d2pos = 0;

		for (size_t i = 0 ; i <= deg ; i++) {
			pos += m_binCoeffV[i] * (uPow[i] * uM1Pow[deg - i]);
		}

		for (size_t i = 0 ; i <= deg ; i++) {
			if (i > 0) {
				dpos += m_binCoeffd1V[i] * (uPow[i - 1] * uM1Pow[deg - i]);
			}
			if (i < deg) {
				dpos -= m_binCoeffd2V[i] * (uPow[i] * uM1Pow[deg - i - 1]);
			}
		}

		for (size_t i = 0 ; i <= deg ; i++) {
			if (i > 1) {
				d2pos += m_binCoeffdd1V[i] * (uPow[i - 2] * uM1Pow[deg - i]);
			}
			if (i > 0 && i < deg) {
				d2pos -= m_binCoeffdd2V[i] * (uPow[i - 1] * uM1Pow[deg - i - 1]);
			}
			if (i < deg - 1) {
				d2pos += m_binCoeffdd3V[i] * (uPow[i] * uM1Pow[deg - i - 2]);
			}
		}
	}
};

