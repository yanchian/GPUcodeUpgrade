#pragma once
#include <iostream>

// Forward declarations
template<typename T> class ScanConvertiblePolygon;

template<int D, typename T>
class Vertex {
};

template<typename T>
class Vertex<2, T> {
public:
	typedef Vector<T, 2> Point;

private:
	Point position_;
 	Point leftNormal_;
	Point rightNormal_;

public:
	Vertex<2, T>(const Point& position, const Point& leftNormal, const Point& rightNormal) :
		position_(position),
		leftNormal_(leftNormal),
		rightNormal_(rightNormal) { }

	Polygon<T> createCharacteristicPolygon(const T range) const;
	ScanConvertiblePolygon<T> createScannableCharacteristicPolygon(const T range) const;

	T characteristicDistance(const Point& p) const;

	Point position() const { return position_; }
	void position(const Point& point) { position_ = point; }

	Point leftNormal() const { return leftNormal_; }
	void leftNormal(const Point& point) { leftNormal_ = point; }

	Point rightNormal() const { return rightNormal_; }
	void rightNormal(const Point& point) { rightNormal_ = point; }
};

