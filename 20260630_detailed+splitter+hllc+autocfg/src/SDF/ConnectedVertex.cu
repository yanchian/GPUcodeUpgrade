#include "ConnectedVertex.hpp"

template<typename T>
void
ConnectedVertex<T>::
init() {
	// get the average normal
	normal_ = 0.0;
	for (typename FaceContainer::const_iterator iter = facesBegin(); iter != facesEnd(); ++iter) {
		normal_ += (*iter)->normal() * (*iter)->angle(this);
	}
	normalise(normal_);

	// establish whether the surface is concave or convex at this point
	int above = 0;
	int below = 0;
	for (typename EdgeContainer::const_iterator iter = edgesBegin(); iter != edgesEnd(); ++iter) {
		const T distance = dot((*iter)->linkedVertex(this)->position() - position(), normal());
		if (distance > 0) {
			above++;
		} else if (distance < 0) {
			below++;
		}
	}
	convex_ = (above == 0);
	concave_ = (below == 0);
}

template<typename T>
Polyhedron<T>
ConnectedVertex<T>::
createCharacteristicPolyhedron(const T range) const {
	if ((convex() && concave()) || faces() < 3) return Polyhedron<T>();

	Polyhedron<T> polyhedron(8);

	// use Gram--Schmidt to generate one vector orthogonal to the normal
	Point orthogonal(0.0);
	if (abs(normal()[0]) < abs(normal()[1])) {
		orthogonal[0] = 1.0;
	} else {
		orthogonal[1] = 1.0;
	}
	orthogonal -= normal() * dot(normal(), orthogonal);
	normalise(orthogonal);

	// generate another using the cross product
	Point orthogonal2 = cross(orthogonal, normal());
	normalise(orthogonal2);

	T minimumDot = std::numeric_limits<T>::max();
	for (typename FaceContainer::const_iterator iter = facesBegin(); iter != facesEnd(); ++iter) {
		minimumDot = std::min(minimumDot, dot((*iter)->normal(), normal()));
	}
	std::cout << "mindot = " << minimumDot << std::endl;

	// TODO: prove this (root 2 tan theta?)
	T tangentLength = sqrt(2.0) * sqrt(1.0 - minimumDot * minimumDot) / minimumDot * range;

	Point average = normal() * range;
	orthogonal *= tangentLength;
	orthogonal2 *= tangentLength;

	Point v[4];
	v[0] = position() + average + orthogonal + orthogonal2;
	v[1] = position() + average + orthogonal - orthogonal2;
	v[2] = position() + average - orthogonal - orthogonal2;
	v[3] = position() + average - orthogonal + orthogonal2;

	for (int i = 0; i < 4; i++) {
		const int ipp = (i + 1) % 4;
		polyhedron.edge(2 * i) = typename Polyhedron<T>::Edge(position(), v[i]);
		polyhedron.edge(2 * i + 1) = typename Polyhedron<T>::Edge(v[i], v[ipp]);
	}

	return polyhedron;
}

template<typename T>
ScanConvertiblePolyhedron<T>
ConnectedVertex<T>::
createScannableCharacteristicPolyhedron(const T range) const {
	return ScanConvertiblePolyhedron<T>(createCharacteristicPolyhedron(range), boost::bind(&ConnectedVertex<T>::characteristicDistance, this, _1), range > 0 ? 1 : -1);
}

template<typename T>
T
ConnectedVertex<T>::
characteristicDistance(ConnectedVertex<T>::Point p) const {
	return abs(p - position());
}

// Defined to sort vertices purely to remove duplicates
template<typename T>
bool operator<(const ConnectedVertex<T>& a, const ConnectedVertex<T>& b) {
	if (a.position()[0] < b.position()[0]) {
		return true;
	}
	if (a.position()[0] > b.position()[0]) {
		return false;
	}
	if (a.position()[1] < b.position()[1]) {
		return true;
	}
	if (a.position()[1] > b.position()[1]) {
		return false;
	}
	if (a.position()[2] < b.position()[2]) {
		return true;
	}
	return false;
}
// Defined to sort vertices purely to remove duplicates
template<typename T>
bool
operator==(const ConnectedVertex<T>& a, const ConnectedVertex<T>& b) {
	for (int i = 0; i < 3; i++) {
		if (a.position()[i] != b.position()[i]) {
			return false;
		}
	}
	return true;
}

template<typename T>
std::ostream&
operator<<(std::ostream& out, const ConnectedVertex<T>& vertex) {
	out << "Vertex {position = " << vertex.position() << ", faces = " << vertex.faces() << ", normal = " << vertex.normal() << ", convex = " << vertex.convex() << ", concave = " << vertex.concave() << "}";
	return out;
}

