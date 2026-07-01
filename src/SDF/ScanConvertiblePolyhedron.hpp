#pragma once
#include "Polyhedron.hpp"
#include "Polygon.hpp"
#include "boost/function.hpp"
#include <algorithm>

// forward declarations
template<typename T> class Polyhedron;

template<typename T>
class ScanConvertiblePolyhedron : Polyhedron<T> {
public:
	typedef typename boost::function<T (const typename Polyhedron<T>::Point)> DistanceFunction;
	typedef Polyhedron<T> Parent;
	
	typedef typename Parent::Container Container;

protected:
	DistanceFunction distance_;
	int sign_;

public:
	ScanConvertiblePolyhedron<T>(const Polyhedron<T>& copy, DistanceFunction distance, const int sign) :
		Polyhedron<T>(copy),
		distance_(distance),
 		sign_(sign) { }

	ScanConvertiblePolyhedron<T>(const typename Polyhedron<T>::Container& points, DistanceFunction distance, const int sign) :
		Polyhedron<T>(points),
		distance_(distance),
		sign_(sign) { }

	using Parent::edgesBegin;
	using Parent::edgesEnd;

	template<typename Grid>
	void scanConvert(Grid& grid) const;
	Polygon<T> cutZSlice(const T z) const;

	int sign() const;
	void sign(const int s);
};


