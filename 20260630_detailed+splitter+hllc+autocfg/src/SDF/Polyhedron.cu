#include "Polyhedron.hpp"

template<typename T>
void
Polyhedron<T>::
insertEdge(const Point& a, const Point& b) {
	edges_.push_back(Edge(a, b));
}


