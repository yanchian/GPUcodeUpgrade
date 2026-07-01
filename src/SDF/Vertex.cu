#include "Vertex.hpp"

template<typename T>
Polygon<T>
Vertex<2, T>::
createCharacteristicPolygon(const T range) const {
	Polygon<T> polygon;

	// make sure this vertex is convex, otherwise it's pointless to return a polygon
	if (range * (leftNormal()[0] * rightNormal()[1] - leftNormal()[1] * rightNormal()[0]) > 0) {
		return polygon;
	}

	polygon.insertPoint(position());
	polygon.insertPoint(position() + range * leftNormal());
	polygon.insertPoint(position() + range * rightNormal());

	return polygon;
}

template<typename T>
T
Vertex<2, T>::
characteristicDistance(const Vertex<2, T>::Point& p) const {
	return abs(p - position());
}

template<typename T>
ScanConvertiblePolygon<T>
Vertex<2, T>::
createScannableCharacteristicPolygon(const T range) const {
	return ScanConvertiblePolygon<T>(createCharacteristicPolygon(range), boost::bind(&Vertex<2, T>::characteristicDistance, this, _1), range > 0 ? 1 : -1);
}

template<typename T>
std::ostream&
operator<<(std::ostream& out, const Vertex<2, T>& vertex) {
	out << "Vertex {position = " << vertex.position() << ", leftNormal = " << vertex.leftNormal() << ", rightNormal = " << vertex.rightNormal() << "}";
	return out;
}

