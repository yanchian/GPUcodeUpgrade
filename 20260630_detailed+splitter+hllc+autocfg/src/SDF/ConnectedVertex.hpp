#pragma once
#include <vector>
#include <iostream>

// Forward declarations
template<typename T> class ConnectedFace;
template<typename T> class ConnectedEdge;
template<typename T> class Polyhedron;
template<typename T> class ScanConvertiblePolyhedron;

/** A vertex which stores adjacency relations with edges and faces.
 *
 * Once initialised, the vertex will have at least one adjacent face and two adjacent edges.
 *
 * \tparam T The numeric type used in representing points.
 */
template<typename T>
class ConnectedVertex {
public:
	typedef Vector<T, 3> Point;
	typedef ConnectedFace<T> Face;
	typedef ConnectedEdge<T> Edge;
	typedef std::vector<Face*> FaceContainer;
	typedef std::vector<Edge*> EdgeContainer;

private:
	Point position_;
	FaceContainer faces_;
	EdgeContainer edges_;

	//@{
	//! Computed variables
	Point normal_;
	bool convex_;
	bool concave_;
	//@}

public:
	//! Trivial constructor
	ConnectedVertex<T>() { }

	//! Once all the adjacency relations are filled, we can populate the computed variables
	void init();

	//@{
	// Accessor functions for the faces associated with this vertex.
	size_t faces() const { return faces_.size(); }
	typename FaceContainer::iterator facesBegin() { return faces_.begin(); }
	typename FaceContainer::iterator facesEnd() { return faces_.end(); }
	typename FaceContainer::const_iterator facesBegin() const { return faces_.begin(); }
	typename FaceContainer::const_iterator facesEnd() const { return faces_.end(); }
	void facesPush(Face* f) { faces_.push_back(f); }
	//@}

	//@{
	// Accessor functions for the edges associated with this vertex.
	size_t edges() const { return edges_.size(); }
	typename EdgeContainer::iterator edgesBegin() { return edges_.begin(); }
	typename EdgeContainer::iterator edgesEnd() { return edges_.end(); }
	typename EdgeContainer::const_iterator edgesBegin() const { return edges_.begin(); }
	typename EdgeContainer::const_iterator edgesEnd() const { return edges_.end(); }
	void edgePush(Edge* e) { edges_.push_back(e); }
	//@}

	//@{
	// Getter and setter for the position of the vertex.
	Point position() const { return position_; }
	void position(const Point& point) { position_ = point; }
	//@}

	//@{
	// Accessors for computed variables
	Point normal() const { return normal_; }
	bool convex() const { return convex_; }
	bool concave() const { return concave_; }
	//@}

	/** The distance to the vertex assuming p is within the characteristic polyhedron.
	 *
	 * \param p A point inside the characteristic polyhedron.
	 */                   
	T characteristicDistance(Point p) const;

	/** Constructs the characteristic polyhedron for this vertex.
	 *
	 *
	 * \param range How far (as a scalar) the polyhedron extends from the boundary.  The sign determines whether the characteristic polyhedron for positive or negative distances is returned.
	 * \return A polyhedron with TODO: complete.
	 */
	Polyhedron<T> createCharacteristicPolyhedron(const T range) const;

	/** Constructs the scannable polyhedron for this vertex.
	 *
	 * \return Same polyhedron as createCharacteristicPolyhedron along with a distance function to the vertex.
	 */
	ScanConvertiblePolyhedron<T> createScannableCharacteristicPolyhedron(const T range) const;
};


