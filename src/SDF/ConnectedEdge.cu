#include "ConnectedEdge.hpp"

template<typename T>
void
ConnectedEdge<T>::
init() {
	if (faces() < 2) return;

	// get a vector in the plane of face 0
	Point b = face(0)->oppositeVertex(this)->position() - vertex(0)->position();

	real s = dot(b, face(1)->normal());

	convex_  = s < 0;
	concave_ = s > 0;

	// construct a normalised tangent vector
	tangent_ = vertex(1)->position() - vertex(0)->position();
	normalise(tangent_);
}

template<typename T>
ConnectedEdge<T>::Vertex*
ConnectedEdge<T>::
linkedVertex(const Vertex* const opposite) const {
	for (int i = 0; i < 2; i++) {
		if (vertices_[i] != opposite) {
			return vertices_[i];
		}
	}
	return NULL;
}

template<typename T>
Polyhedron<T>
ConnectedEdge<T>::
createCharacteristicPolyhedron(T range) const {
	if (vertices() != 2 || faces() != 2) return Polyhedron<T>();
	Polyhedron<T> polyhedron(9);

	// get the angle between the adjacent normals; this is used to scale the
	// range so that the polyhedron encompasses the entire domain
	T cosTheta = dot(face(0)->normal(), face(1)->normal());
	range /= cosTheta;
	int n = 0;

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			polyhedron.edge(n++) = typename Polyhedron<T>::Edge(vertex(j)->position(), vertex(j)->position() + range * face(i)->normal());
		}
		polyhedron.edge(n++) = typename Polyhedron<T>::Edge(vertex(0)->position() + range * face(i)->normal(), vertex(1)->position() + range * face(i)->normal());
		polyhedron.edge(n++) = typename Polyhedron<T>::Edge(vertex(i)->position() + range * face(0)->normal(), vertex(i)->position() + range * face(1)->normal());
	}
	polyhedron.edge(n++) = typename Polyhedron<T>::Edge(vertex(0)->position(), vertex(1)->position());

	return polyhedron;
}

template<typename T>
ScanConvertiblePolyhedron<T>
ConnectedEdge<T>::
createScannableCharacteristicPolyhedron(const T range) const {
	return ScanConvertiblePolyhedron<T>(createCharacteristicPolyhedron(range), boost::bind(&ConnectedEdge<T>::characteristicDistance, this, _1), range > 0 ? 1 : -1);
}

template<typename T>
T
ConnectedEdge<T>::
characteristicDistance(ConnectedEdge<T>::Point p) const {
	return abs(cross(p - vertex(0)->position(), tangent()));
}

// Defined to sort edges purely to remove duplicates
template<typename T>
bool
operator<(const ConnectedEdge<T>& a, const ConnectedEdge<T>& b) {
	const typename ConnectedEdge<T>::Vertex* va1;
	const typename ConnectedEdge<T>::Vertex* va2;
	const typename ConnectedEdge<T>::Vertex* vb1;
	const typename ConnectedEdge<T>::Vertex* vb2;

	if ((a.vertex(0)) < (a.vertex(1))) {
		va1 = a.vertex(0);
		va2 = a.vertex(1);
	} else {
		va1 = a.vertex(1);
		va2 = a.vertex(0);
	}
	if ((b.vertex(0)) < (b.vertex(1))) {
		vb1 = b.vertex(0);
		vb2 = b.vertex(1);
	} else {
		vb1 = b.vertex(1);
		vb2 = b.vertex(0);
	}

	if (va1 < vb1) return true;
	if (va1 > vb1) return false;
	if (va2 < vb2) return true;
	return false;
}

template<typename T>
bool
operator==(const ConnectedEdge<T>& a, const ConnectedEdge<T>& b) {
	return (a.vertex(0) == b.vertex(0) && a.vertex(1) == b.vertex(1)) || (a.vertex(0) == b.vertex(1) && a.vertex(1) == b.vertex(0));
}

template<typename T>
std::ostream&
operator<<(std::ostream& out, const ConnectedEdge<T>& edge) {
	out << "Edge {start = " << edge.vertex(0) << " " << edge.vertex(0)->position() << ", end = " << edge.vertex(1) << " " << edge.vertex(1)->position() << ", faces = " << edge.faces() << ", convex = " << edge.convex() << ", concave = " << edge.concave() << "}";
	return out;
}

