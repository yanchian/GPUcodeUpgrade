#pragma once
#include <iostream>
#include "Polygon.hpp"
#include "ScanConvertiblePolygon.hpp"

// Forward declarations
template<typename T> class Polyhedron;

template<int D, typename T>
class Face {
};

template<typename T>
class Face<2, T> {
public:
	typedef Vector<T, 2> Point;

private:
	Point start_;
	Point end_;
	Point normal_;
	Point tangent_;

public:
	Face<2, T>(const Point& start, const Point& end, const Point& normal) :
		start_(start),
		end_(end),
		normal_(normal) {
		tangent_ = end - start;
		tangent_ /= abs(tangent_);
	}

	//! Build a crude characteristic polygon by casting out by !range
	Polygon<T> createCharacteristicPolygon(const T range) const;
	ScanConvertiblePolygon<T> createScannableCharacteristicPolygon(const T range) const;

	//! Return the distance to the edge, assuming we're inside the characteristic polygon
	T characteristicDistance(Point p) const;

	Point start() const { return start_; }
	void start(const Point& point) { start_ = point; }

	Point end() const { return end_; }
	void end(const Point& point) { end_ = point; }

	Point normal() const { return normal_; }
	void normal(const Point& point) { normal_ = point; }

	Point tangent() const { return tangent_; }
	void tangent(const Point& point) { tangent_ = point; }
};

