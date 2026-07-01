#pragma once
#include "Polygon.hpp"
#include "boost/function.hpp"
#include <algorithm>

template<typename T>
class ScanConvertiblePolygon : Polygon<T> {
public:
	typedef typename boost::function<T (const typename Polygon<T>::Point)> DistanceFunction;
	typedef Polygon<T> Parent;
protected:
	DistanceFunction distance_;
	int sign_;

public:
	ScanConvertiblePolygon<T>(const Polygon<T>& copy, DistanceFunction distance, const int sign) :
		Polygon<T>(copy),
		distance_(distance),
 		sign_(sign) { }

	ScanConvertiblePolygon<T>(const typename Polygon<T>::Container& points, DistanceFunction distance, const int sign) :
		Polygon<T>(points),
		distance_(distance),
		sign_(sign) { }

	using Parent::points;
	using Parent::pointsBegin;
	using Parent::pointsEnd;

	template<typename Grid, typename Writer>
	void scanConvert(Grid& grid, Writer& writer) const;
	template<typename Grid>
	void scanConvert(Grid& grid) const;

	std::vector<T> getXIntersections(const T y) const;

	int sign() const;
	void sign(const int s);
};


