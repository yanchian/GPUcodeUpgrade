#pragma once

template<typename T>
Polygon<T>
Face<2, T>::
createCharacteristicPolygon(const T range) const {
	Polygon<T> polygon;

	Point edge = range * normal();

	polygon.insertPoint(start());
	polygon.insertPoint(start() + edge);
	polygon.insertPoint(end() + edge);
	polygon.insertPoint(end());

	return polygon;
}

template<typename T>
ScanConvertiblePolygon<T>
Face<2, T>::
createScannableCharacteristicPolygon(const T range) const {
	return ScanConvertiblePolygon<T>(createCharacteristicPolygon(range), boost::bind(&Face<2, T>::characteristicDistance, this, _1), range > 0 ? 1 : -1);
}

template<typename T>
T
Face<2, T>::
characteristicDistance(Face<2, T>::Point p) const {
	p -= start();
	return std::abs(p[0] * tangent()[1] - p[1] * tangent()[0]);
}

template<typename T>
std::ostream&
operator<<(std::ostream& out, const Face<2, T>& face) {
	out << "Face {start = " << face.start() << ", end = " << face.end() << ", normal = " << face.normal() << ", tangent = " << face.tangent() << "}";
	return out;
}

