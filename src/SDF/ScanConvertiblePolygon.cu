template<typename Grid, typename T>
class Grid2DWriter {
private:
	Grid& grid;
	typedef typename boost::function<T (const typename Polygon<T>::Point)> DistanceFunction;
	DistanceFunction func;

public:
	Grid2DWriter(Grid& grid_, DistanceFunction func_) :
		grid(grid_),
 		func(func_)	{ }

	typename Grid::Type&
	operator()(const int i, const int j) {
		return grid(i, j, 0);
	}

	bool
	exists(const int i, const int j) {
		return grid.exists(i, j);
	}

	T
	distance(const int i, const int j) {
		return func(typename Polygon<T>::Point(grid.x(i), grid.y(j)));
	}
};

template<typename T>
template<typename Grid>
void
ScanConvertiblePolygon<T>::
scanConvert(Grid& grid) const {
	Grid2DWriter<Grid, T> writer(grid, distance_);
	scanConvert(grid, writer);
}

template<typename T>
template<typename Grid, typename Writer>
void
ScanConvertiblePolygon<T>::
scanConvert(Grid& grid, Writer& writer) const {
	// get the extent of the polygon in y
	T minimumY = std::numeric_limits<T>::infinity();
	T maximumY = -std::numeric_limits<T>::infinity();
	for (typename Polygon<T>::ConstIterator iter = Polygon<T>::pointsBegin(); iter != Polygon<T>::pointsEnd(); ++iter) {
		T y = (*iter)[1];
		if (y < minimumY) {
			minimumY = y;
		}
		if (y > maximumY) {
			maximumY = y;
		}
	}

	// clip to the size of the grid
	minimumY = std::max(minimumY, grid.y(-grid.ghostCells()));
	maximumY = std::min(maximumY, grid.y(grid.activeNy() + grid.ghostCells() - 1));

	int minimumJ = std::max<int>(std::ceil(grid.j(minimumY)), -grid.ghostCells());
	int maximumJ = std::min<int>(std::floor(grid.j(maximumY)), grid.activeNy() + grid.ghostCells() - 1);

	for (int j = minimumJ; j <= maximumJ; j++) {
		std::vector<T> intersections = getXIntersections(grid.y(j));
		//std::cout << intersections.size() << " intersections at " << grid.y(j) << std::endl;
		std::sort(intersections.begin(), intersections.end());

		if (intersections.size() % 2 != 0) {
			std::cout << "oops: " << grid.y(j) << " ";
			for (int i = 0; i < intersections.size(); i++) {
				std::cout << intersections[i] << " ";
			}
			std::cout << std::endl << *this << std::endl;
			std::cout << "OH NO, OH GOD NO." << std::endl;
			return;
		}
		if (intersections.size() > 0) {
			for (int k = 0; k < intersections.size(); k += 2) {
				const int minimumI = std::max<int>(std::ceil(grid.i(intersections[k])), -grid.ghostCells());
				const int maximumI = std::min<int>(std::floor(grid.i(intersections[k + 1])), grid.activeNx() + grid.ghostCells() - 1);

				for (int i = minimumI; i <= maximumI; i++) {
					const real d = writer.distance(i, j);
					if (d < abs(writer(i, j))) {
						writer(i, j) = sign() * d;
					}
				}
			}
		}
	}
}

template<typename T>
std::vector<T>
ScanConvertiblePolygon<T>::
getXIntersections(const T y) const {
	std::vector<T> intersections;
	for (int i = 0; i < points(); i++) {
		const typename Polygon<T>::Point start = Polygon<T>::point(i);
		const typename Polygon<T>::Point end   = Polygon<T>::point(i + 1);

		if (
			((start[1] < y) != (end[1] <= y) && start[1] != y && end[1] != y && end[1] != start[1]) ||
			(end[1] == y && (end[1] < start[1]) != (end[1] <= Polygon<T>::point(i + 2)[1]))
			) {
			intersections.push_back((start[0] * (end[1] - y) + end[0] * (y - start[1])) / (end[1] - start[1]));
		}
	}
	return intersections;
}

template<typename T>
int
ScanConvertiblePolygon<T>::
sign() const {
	return sign_;
}

template<typename T>
void
ScanConvertiblePolygon<T>::
sign(const int s) {
	sign_ = s;
}


