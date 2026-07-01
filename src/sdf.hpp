#include "StridedArray.hpp"

StridedCell<real, 2> typedef LS;

class SDFUnaryOperation {
public:
	virtual ~SDFUnaryOperation() { }
	virtual void operator()(const LS& a, LS r) const = 0;
};

class SDFInvert : public SDFUnaryOperation {
public:
	virtual void operator()(const LS& a, LS r) const {
		r[0] = -a[0];
		r[1] = 1.0 - a[1];
	}
};

class SDFBinaryOperation {
public:
	virtual ~SDFBinaryOperation() { }
	virtual void operator()(const LS& a, const LS& b, LS r) const = 0;
};

class SDFUnion : public SDFBinaryOperation {
public:
	virtual void operator()(const LS& a, const LS& b, LS r) const {
		if (a[0] > b[0]) {
			r[0] = a[0];
			r[1] = a[1];
		} else {
			r[0] = b[0];
			r[1] = b[1];
		}
	}
};

class SDFIntersection : public SDFBinaryOperation {
public:
	virtual void operator()(const LS& a, const LS& b, LS r) const {
		if (a[0] < b[0]) {
			r[0] = a[0];
			r[1] = a[1];
		} else {
			r[0] = b[0];
			r[1] = b[1];
		}
	}
};

class Shape {
public:
	virtual ~Shape() { }
	virtual void distance(Vec p, LS r, Vec d) const = 0;
};

class Rectangle : public Shape {
private:
	Vec corners_[2];

public:
	Rectangle(Vec corners[2]) {
		for (int i = 0; i < 2; i++) corners_[i] = corners[i];
	}

	void distance(Vec p, LS r, Vec d) const {
		if (p[0] > corners_[0][0] && p[1] > corners_[0][1] && p[0] < corners_[1][0] && p[1] < corners_[1][1]) {
			// inside rectangle
			r[0] = min(corners_[1][0] - p[0], min(p[0] - corners_[0][0], min(corners_[1][1] - p[1], p[1] - corners_[0][1])));
		} else {
			// outside rectangle
			if (p[0] < corners_[0][0]) {
				if (p[1] > corners_[1][1]) {
					r[0] = -sqrt(pow(p[0] - corners_[0][0], 2) + pow(p[1] - corners_[1][1], 2));
				} else if (p[1] < corners_[0][1]) {
					r[0] = -sqrt(pow(p[0] - corners_[0][0], 2) + pow(p[1] - corners_[0][1], 2));
				} else {
					r[0] = -(corners_[0][0] - p[0]);
				}
			} else if (p[0] > corners_[1][0]) {
				if (p[1] > corners_[1][1]) {
					r[0] = -sqrt(pow(p[0] - corners_[1][0], 2) + pow(p[1] - corners_[1][1], 2));
				} else if (p[1] < corners_[0][1]) {
					r[0] = -sqrt(pow(p[0] - corners_[1][0], 2) + pow(p[1] - corners_[0][1], 2));
				} else {
					r[0] = -(p[0] - corners_[1][0]);
				}
			} else if (p[1] > corners_[1][1]) {
				r[0] = -(p[1] - corners_[1][1]);
			} else {
				r[0] = -(corners_[0][1] - p[1]);
			}
		}
	}
};

class Circle : public Shape {
private:
	Vec centre_;
	real radius_;

public:
	Circle(Vec centre, real radius) : centre_(centre), radius_(radius)
	{ }

	void distance(Vec p, LS r, Vec d) const {
		r[0] = radius_ - norm(p - centre_);
		const int N = 32;
		if (r[0] < 0.5 * norm(d)) {
			r[1] = 0.0;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					const real x = p[0] - centre_[0] + (i - N / 2) * d[0] / N;
					const real y = p[1] - centre_[1] + (j - N / 2) * d[1] / N;
					if (x * x + y * y < radius_ * radius_) {
						r[1] += 1.0 / (N * N);
					}
				}
			}
		} else {
			r[1] = r[0] > 0;
		}
	}
};

__host__ __device__ Vec rotateCoordinates(const Vec p, const Vec o, const real theta) {
	Vec r;
	r[0] = (p[0] - o[0]) * cos(theta) - (p[1] - o[1]) * sin(theta) + o[0];
	r[1] = (p[0] - o[0]) * sin(theta) + (p[1] - o[1]) * cos(theta) + o[1];
	return r;
}

__host__ __device__ real pointLineDistance(real x, real y, real ax, real ay, real bx, real by) {
	x -= ax;
	y -= ay;
	bx -= ax;
	by -= ay;
	const real pb = x * bx + y * by;
	const real bb = bx * bx + by * by;
	if (pb > bb) {
		x -= bx;
		y -= by;
		return -sqrt(x * x + y * y);
	} else if (pb < 0) {
		return -sqrt(x * x + y * y);
	} else {
		return (x * by - y * bx) / sqrt(bb);
	}
}

