#include "ScanConvertiblePolyhedron.hpp"

template<typename Grid, typename T>
class Grid3DWriter {
private:
	Grid& grid;
	const int k;
	typedef typename boost::function<T (const typename Polyhedron<T>::Point)> DistanceFunction;
	DistanceFunction func;

public:
	Grid3DWriter(Grid& grid_, DistanceFunction func_, const int k_) :
		grid(grid_),
		k(k_),
 		func(func_)	{ }

	typename Grid::Type&
	operator()(const int i, const int j) {
		return grid(i, j, k, 0);
	}

	bool
	exists(const int i, const int j) {
		return grid.active(i, j, k);
	}

	T
	distance(const int i, const int j) {
		return func(typename Polyhedron<T>::Point(grid.x(i), grid.y(j), grid.z(k)));
	}
};

template<typename T>
template<typename Grid>
void
ScanConvertiblePolyhedron<T>::
scanConvert(Grid& grid) const {
	// get the extent of the polygon in z
	T minimumZ = std::numeric_limits<T>::infinity();
	T maximumZ = -std::numeric_limits<T>::infinity();
	for (typename Container::const_iterator iter = edgesBegin(); iter != edgesEnd(); ++iter) {
		minimumZ = std::min((*iter).start()[2], minimumZ);
		maximumZ = std::max((*iter).start()[2], maximumZ);
		minimumZ = std::min((*iter).end()[2], minimumZ);
		maximumZ = std::max((*iter).end()[2], maximumZ);
	}

	// clip to the size of the grid
	minimumZ = std::max(minimumZ, grid.z(0));
	maximumZ = std::min(maximumZ, grid.z(grid.activeNz() - 1));

	int minimumK = grid.k(minimumZ) - 0.5;
	int maximumK = grid.k(maximumZ) + 1.0;

	for (int k = minimumK; k < maximumK; k++) {
		Polygon<T> polygon = cutZSlice(grid.z(k));
		polygon.orderVertices();
		//std::cout << "at z=" << k << " (" << sign_ << ") " << grid.z(k) << " " << polygon << std::endl;
		Grid3DWriter<Grid, T> writer(grid, distance_, k);
		ScanConvertiblePolygon<T>(polygon, distance_, sign_).scanConvert(grid, writer);
	}
}

template<typename T>
Polygon<T>
ScanConvertiblePolyhedron<T>::
cutZSlice(const T z) const {
	Polygon<T> polygon;
	for (typename Container::const_iterator iter = edgesBegin(); iter != edgesEnd(); ++iter) {
		if (((*iter).start()[2] < z) ^ ((*iter).end()[2] <= z)) {
			typename Polygon<T>::Point p;
			const real t = ((*iter).end()[2] - z) / ((*iter).end()[2] - (*iter).start()[2]);
			for (int i = 0; i < 2; i++) {
				p[i] = (*iter).start()[i] * t + (*iter).end()[i] * (1.0 - t);
			}
			polygon.insertPoint(p);
		}
	}
	return polygon;
}

template<typename T>
int
ScanConvertiblePolyhedron<T>::
sign() const {
	return sign_;
}

template<typename T>
void
ScanConvertiblePolyhedron<T>::
sign(const int s) {
	sign_ = s;
}


