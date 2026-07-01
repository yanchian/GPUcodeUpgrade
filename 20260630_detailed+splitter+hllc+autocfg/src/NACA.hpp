/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#pragma once
#include "SDF/Polygon.hpp"

class NACAAirfoil {
private:
	const real c_, t_, m_, p_;

	const real A, B, C, D, E;
public:
	NACAAirfoil(const real c, const real t, const real m, const real p) : c_(c), t_(t), m_(m), p_(p), A(0.2969), B(-0.1260), C(-0.3516), D(0.2843), E(-0.1036) {
	}

	real y_t(const real x) const {
		const real xc = x / c_;
		return 5.0 * t_ * c_ * (A * sqrt(xc) + B * xc + C * xc * xc + D * xc * xc * xc + E * xc * xc * xc * xc);
	}
	real y_c(const real x) const {
		if (0.0 <= x && x <= p_ * c_) {
			return (m_ * x) / (p_ * p_) * (2.0 * p_ - x / c_);
		} else {
			return m_ * (c_ - x) / pow(1.0 - p_, 2.0) * (1.0 + x / c_ - 2.0 * p_);
		}
	}
	real dy_c(const real x) const {
		if (0.0 <= x && x <= p_ * c_) {
			return (m_ * (2.0 * p_ - x / c_))/(p_ * p_) - (m_ * x)/(c_ * p_ * p_);
		} else {
			return (m_ * (c_ - x))/(c_ * pow(1.0 - p_, 2.0)) - (m_ * (1.0 + x / c_ - 2.0 * p_))/pow(1.0 - p_, 2.0);
		}
	}

	real xTop(const real l) const {
		return l - y_t(l) * sin(atan(dy_c(l)));
	}
	real yTop(const real l) const {
		return y_c(l) + y_t(l) * cos(atan(dy_c(l)));
	}
	real xBottom(const real l) const {
		return l + y_t(l) * sin(atan(dy_c(l)));
	}
	real yBottom(const real l) const {
		return y_c(l) - y_t(l) * cos(atan(dy_c(l)));
	}
	real c() const { return c_; }
	real m() const { return m_; }
	real p() const { return p_; }

	Polygon<real> getApproximation(const int n) const {
		Polygon<real> polygon(0);

		const real delta = 2.0 * c_ / n;
		for (int i = 0; i < n; i++) {
			real l = i * delta;
			if (l <= c_) {
				polygon.insertPoint(Polygon<real>::Point(xBottom(l), yBottom(l)));
			} else {
				l = 2.0 * c_ - l;
				polygon.insertPoint(Polygon<real>::Point(xTop(l), yTop(l)));
			}
		}
		return polygon;
	}
};

