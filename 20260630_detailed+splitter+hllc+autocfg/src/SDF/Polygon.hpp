#pragma once
#include "../Vector.hpp"
#include <vector>
#include <cstddef>

template<typename T>
class Polygon {
public:
	typedef Vector<T, 2> Point;
	typedef std::vector<Point> Container;
	typedef typename Container::iterator Iterator;
	typedef typename Container::const_iterator ConstIterator;

protected:
	Container points_;

public:
	Polygon<T>() { }

	Polygon<T>(const Container& points) :
		points_(points) { }

	Polygon<T>(const size_t size) :
		points_(size) { }

	Iterator pointsBegin() { return points_.begin(); }
	Iterator pointsEnd() { return points_.end(); }

	ConstIterator pointsBegin() const { return points_.begin(); }
	ConstIterator pointsEnd() const { return points_.end(); }

	int points() const { return points_.size(); }
	Point point(int i) const;
	Point& point(int i);
	void insertPoint(const Point& p) { points_.push_back(p); }

	void orderVertices();
};

template<typename T>
Polygon<T>::Point
Polygon<T>::
point(int i) const {
	i = i % (int)points_.size();
	if (i < 0) i += points_.size();
	return points_[i];
}

template<typename T>
Polygon<T>::Point&
Polygon<T>::
point(int i) {
	i = i % (int)points_.size();
	if (i < 0) i += points_.size();
	return points_[i];
}

template<typename T>
struct PointWithAngle {
	Vector<T, 2> point;
	T angle;
};

template<typename T>
bool
operator<(const PointWithAngle<T>& a, const PointWithAngle<T>& b) {
	return a.angle < b.angle;
}

template<typename T>
void
Polygon<T>::
orderVertices() {
	if (points() <= 3) return;
	// find the point with smallest x coordinate and move it to the first place in the array
	int first = 0;
	for (int i = 1; i < points(); i++) {
		if (point(i)[0] < point(first)[0]) {
			first = i;
		}
	}
	if (first != 0) {
		std::swap(point(0), point(first));
	}

	std::vector<PointWithAngle<T> > angles(points() - 1);
	for (int i = 1; i < points(); i++) {
		Point rel = point(i) - point(0);
		PointWithAngle<T> angle = {points_[i], rel[1] / abs(rel)};
		// if the point is coincident, make sure it comes first in the array
		if (abs(rel) == 0.0) {
			angle.angle = -1.0;
		}
		angles[i - 1] = angle;
	}
	std::sort(angles.begin(), angles.end());
	for (int i = 1; i < points(); i++) {
		points_[i] = angles[i - 1].point;
	}
}

template<typename T>
std::ostream&
operator<<(std::ostream& out, const Polygon<T>& polygon) {
	out << "Polygon {points: [";
	for (typename Polygon<T>::Container::const_iterator iter = polygon.pointsBegin(); iter != polygon.pointsEnd(); ++iter) {
		out << *iter << ", ";
	}
	out << "]}";
	return out;
}

