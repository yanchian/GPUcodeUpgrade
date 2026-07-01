#include "ConnectedFace.hpp"

template<typename T>
void
ConnectedFace<T>::
init() {
	area_ = 0.5 * abs(cross(vertex(1)->position() - vertex(0)->position(), vertex(2)->position() - vertex(0)->position()));
	normalise(normal_);
}

template<typename T>
T
ConnectedFace<T>::
angle(const ConnectedFace<T>::Vertex * const a) const {
	// find the other two vertices
	const Vertex * b = NULL;
	const Vertex * c = NULL;
	for (int i = 0; i < 3; i++) {
		if (a != vertices_[i]) {
			if (b == NULL) {
				b = vertices_[i];
			} else {
				c = vertices_[i];
			}
		}
	}

	if (a == NULL || b == NULL) {
		return 0.0;
	} else {
		return vectorAngle(b->position() - a->position(), c->position() - a->position());
	}
}

template<typename T>
T
ConnectedFace<T>::
oppositeMidpoint(const ConnectedFace<T>::Vertex * const a) const {
	// find the other two vertices
	const Vertex * b = NULL;
	const Vertex * c = NULL;
	for (int i = 0; i < 3; i++) {
		if (a != vertices_[i]) {
			if (b == NULL) {
				b = vertices_[i];
			} else {
				c = vertices_[i];
			}
		}
	}

	return 0.5 * (b->position() + c->position());
}

template<typename T>
const ConnectedFace<T>::Vertex*
ConnectedFace<T>::
oppositeVertex(const ConnectedFace<T>::Edge * const edge) const {
	for (int i = 0; i < 3; i++) {
		if (vertices_[i] != edge->vertex(0) && vertices_[i] != edge->vertex(1)) {
			return vertices_[i];
		}
	}
	return NULL;
}
template<typename T>
Polyhedron<T>
ConnectedFace<T>::
createCharacteristicPolyhedron(const T range) const {
	Polyhedron<T> polyhedron(9);
	int n = 0;

	Point edge = range * normal();

	for (int i = 0; i < 3; i++) {
		const int ip = (i + 1) % 3;
		polyhedron.edge(n++) = typename Polyhedron<T>::Edge(vertex(i)->position(), vertex(ip)->position());
		polyhedron.edge(n++) = typename Polyhedron<T>::Edge(vertex(i)->position() + edge, vertex(ip)->position() + edge);
		polyhedron.edge(n++) = typename Polyhedron<T>::Edge(vertex(i)->position(), vertex(i)->position() + edge);
	}

	return polyhedron;
}

template<typename T>
ScanConvertiblePolyhedron<T>
ConnectedFace<T>::
createScannableCharacteristicPolyhedron(const T range) const {
	return ScanConvertiblePolyhedron<T>(createCharacteristicPolyhedron(range), boost::bind(&ConnectedFace<T>::characteristicDistance, this, _1), range > 0 ? 1 : -1);
}

template<typename T>
T
ConnectedFace<T>::
characteristicDistance(ConnectedFace<T>::Point p) const {
	return abs(dot(normal(), p - vertex(0)->position()));
}

template<typename T>
std::ostream&
operator<<(std::ostream& out, const ConnectedFace<T>& face) {
	out << "Face {normal = " << face.normal() << ", vertices = [";
	for (typename ConnectedFace<T>::VertexContainer::const_iterator iter = face.verticesBegin(); iter != face.verticesEnd(); ++iter) {
		out << (*iter)->position() << ", ";
	}
	out << "]}";
	return out;
}
