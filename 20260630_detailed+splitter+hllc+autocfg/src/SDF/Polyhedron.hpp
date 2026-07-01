#pragma once
#include <vector>
#include <cstddef>
#include "Edge.hpp"

// Forward declarations

template<typename T>
class Polyhedron {
public:
	typedef Vector<T, 3> Point;
	typedef Edge<3, T> Edge;
	typedef std::vector<Edge> Container;

protected:
	Container edges_;

public:
	Polyhedron<T>() :
		edges_(0) { }

	Polyhedron<T>(const Container& edges) :
		edges_(edges) { }

	Polyhedron<T>(const size_t size) :
		edges_(size) { }

	size_t edges() const { return edges_.size(); }
	typename Container::iterator edgesBegin() { return edges_.begin(); }
	typename Container::iterator edgesEnd() { return edges_.end(); }
	typename Container::const_iterator edgesBegin() const { return edges_.begin(); }
	typename Container::const_iterator edgesEnd() const { return edges_.end(); }
	//Edge  edge(const int i) const { return edges_[i]; }
	Edge& edge(const int i) { return edges_[i]; }
	void insertEdge(const Edge& e) { edges_.push_back(e); }
	void insertEdge(const Point& a, const Point& b);
};


