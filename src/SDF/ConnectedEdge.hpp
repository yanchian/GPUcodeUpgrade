#pragma once
#include <vector>
#include <iostream>

// Forward declarations
template<typename T> class ConnectedVertex;
template<typename T> class ConnectedFace;
template<typename T> class Polyhedron;
template<typename T> class ScanConvertiblePolyhedron;

/** A vertex which stores adjacency relations with faces and vertices.
 *
 * Once initialised, the edge will have exactly two associated vertices and at least one face.
 *
 * \tparam T The numeric type used in representing points.
 */
template<typename T>
class ConnectedEdge {
public:
	typedef Vector<T, 3> Point;
	typedef ConnectedVertex<T> Vertex;
	typedef ConnectedFace<T> Face;

	typedef std::vector<Vertex*> VertexContainer;
	typedef std::vector<Face*> FaceContainer;

private:
	VertexContainer vertices_;
	FaceContainer faces_;

	//@{
	// Computed variables
	Point tangent_;
	bool convex_;
	bool concave_;
	//@}

public:
	ConnectedEdge<T>() : vertices_(2) { };

	//! Once all the adjacency relations are filled, we can populate the computed variables
	void init();

	//@{
	// Accessor functions for the vertices associated with this edge.
	//
	// Once initialised, it may be assumed that the edge has exactly two vertices.
	size_t vertices() const { return vertices_.size(); }
	typename VertexContainer::iterator verticesBegin() { return vertices_.begin(); }
	typename VertexContainer::iterator verticesEnd() { return vertices_.end(); }
	typename VertexContainer::const_iterator verticesBegin() const { return vertices_.begin(); }
	typename VertexContainer::const_iterator verticesEnd() const { return vertices_.end(); }
	const Vertex* vertex(const int i) const { return vertices_[i]; }
	Vertex*& vertex(const int i) { return vertices_[i]; }
	//@}

	//@{
	// Accessor functions for faces joined to this edge
	//
	// Once initialised, it may not be assumed that the edge has exactly two adjacent faces.  If the mesh is not well-formed (closed), or is self-intersecting, it may have either more or fewer.
	size_t faces() const { return faces_.size(); }
	typename FaceContainer::iterator facesBegin() { return faces_.begin(); }
	typename FaceContainer::iterator facesEnd() { return faces_.end(); }
	typename FaceContainer::const_iterator facesBegin() const { return faces_.cbegin(); }
	typename FaceContainer::const_iterator facesEnd() const { return faces_.cend(); }
	const Face* face(const int i) const { return faces_[i]; }
	Face*& face(const int i) { return faces_[i]; }
	void facesPushBack(Face* f) { faces_.push_back(f); }

	//! Retrieves a pointer to the edge vertex other than the one given.
	Vertex* linkedVertex(const Vertex* const opposite) const;

	//@{
	// Accessors for computed variables
	bool convex() const { return convex_; }
	bool concave() const { return concave_; }
	Point tangent() const { return tangent_; }
	//@}

	/** The distance to the edge assuming p is within the characteristic polyhedron.
	 *
	 * If p is not inside the characteristic polyhedron, the distance may be an underestimate.
	 *
	 * \param p A point inside the characteristic polyhedron.
	 */                   
	T characteristicDistance(Point p) const;

	/** Constructs the characteristic polyhedron for this edge.
	 *
	 *
	 * \param range How far (as a scalar) the polyhedron extends from the boundary.  The sign determines whether the characteristic polyhedron for positive or negative distances is returned.
	 * \return A polyhedron with nine edges if the edge is locally closed (joins two faces); an empty polyhedron otherwise.
	 */
	Polyhedron<T> createCharacteristicPolyhedron(T range) const;

	/** Constructs the scannable polyhedron for this edge.
	 *
	 * \return Same polyhedron as createCharacteristicPolyhedron along with a distance function to the edge.
	 */
	ScanConvertiblePolyhedron<T> createScannableCharacteristicPolyhedron(const T range) const;
};


