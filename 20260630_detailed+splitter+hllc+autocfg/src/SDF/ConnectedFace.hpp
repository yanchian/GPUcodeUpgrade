#pragma once
#include <vector>
#include <iostream>

// Forward declarations
template<typename T> class ConnectedEdge;
template<typename T> class ConnectedVertex;
template<typename T> class Polyhedron;
template<typename T> class ScanConvertiblePolyhedron;

/** A face which stores adjacency relations with edges.
 *
 * Once initialised, the face is guaranteed to have three vertices and three edges.
 *
 * \tparam T The numeric type used in representing points.
 */
template<typename T>
class ConnectedFace {
public:
	typedef Vector<T, 3> Point;
	typedef ConnectedEdge<T> Edge;
	typedef ConnectedVertex<T> Vertex;

	typedef std::vector<Edge*> EdgeContainer;
	typedef std::vector<Vertex*> VertexContainer;

private:
	EdgeContainer edges_;
	VertexContainer vertices_;
	Point normal_;

	//@{
	// Computed variables
	T area_;
	//@}

public:
	ConnectedFace<T>() : vertices_(3), edges_(3) { }

	//! Once all the adjacency relations are filled, we can populate the computed variables
	void init();

	//@{
	// Accessor functions for the vertices associated with this face.
	//
	// Once initialised, it may be assumed that the edge has exactly three vertices
	size_t vertices() const { return vertices_.size(); }
	typename VertexContainer::iterator verticesBegin() { return vertices_.begin(); }
	typename VertexContainer::iterator verticesEnd() { return vertices_.end(); }
	typename VertexContainer::const_iterator verticesBegin() const { return vertices_.begin(); }
	typename VertexContainer::const_iterator verticesEnd() const { return vertices_.end(); }
	Vertex const * vertex(const int i) const { return vertices_[i]; }
	Vertex*& vertex(const int i) { return vertices_[i]; }
	//@}

	//@{
	// Accessor functions for the vertices associated with this face.
	//
	// Once initialised, it may be assumed that the edge has exactly three edges.
	size_t edges() const { return edges_.size(); }
	typename EdgeContainer::iterator edgesBegin() { return edges_.begin(); }
	typename EdgeContainer::iterator edgesEnd() { return edges_.end(); }
	typename EdgeContainer::const_iterator edgesBegin() const { return edges_.begin(); }
	typename EdgeContainer::const_iterator edgesEnd() const { return edges_.end(); }
	Edge const * edge(const int i) const { return edges_[i]; }
	Edge*& edge(const int i) { return edges_[i]; }
	//@}

	//@{
	// Accessors for computed variables
	Point normal() const { return normal_; }
	void normal(const Point& p) { normal_ = p; }

	T area() const { return area_; }
	//@}

	/** Get the angle around a vertex.
	 *
	 * \param a Pointer to the vertex around which we want the angle.
	 */
	T angle(const Vertex * const a) const;

	/** Get the midpoint of the edge opposite a.
	 */
	T oppositeMidpoint(const Vertex * const a) const;

	/** Find the vertex opposite an edge.
	 *
	 * \param edge An edge which is associated with this face.
	 * \return A pointer to the vertex on this edge which is not on the given edge.
	 */
	const Vertex* oppositeVertex(const Edge * const edge) const;

	/** The distance to the face assuming p is within the characteristic polyhedron.
	 *
	 * If p is not inside the polyhedron, the distance may be an underestimate.
	 *
	 * \param p A point inside the characteristic polyhedron.
	 */                   
	T characteristicDistance(Point p) const;

	/** Constructs the characteristic polyhedron for this face.
	 *
	 *
	 * \param range How far (as a scalar) the polyhedron extends from the boundary.  The sign determines whether the characteristic polyhedron for positive or negative distances is returned.
	 * \return A polyhedron with nine edges if the face is locally closed (joins two faces); an empty polyhedron otherwise.
	 */
	Polyhedron<T> createCharacteristicPolyhedron(const T range) const;

	/** Constructs the scannable polyhedron for this face.
	 *
	 * \return Same polyhedron as createCharacteristicPolyhedron along with a distance function to the face.
	 */
	ScanConvertiblePolyhedron<T> createScannableCharacteristicPolyhedron(const T range) const;
};


