#pragma once
#include "../Vector.hpp"
#include "Polygon.hpp"
#include "Face.hpp"
#include "Vertex.hpp"
#include "ConnectedFace.hpp"
#include "ConnectedEdge.hpp"
#include "ConnectedVertex.hpp"
#include <vector>
#include <fstream>
#include <boost/cstdint.hpp>

// Forward declarations

/** \addtogroup AvailableGeom Geometries
    <H2>%STL format</H2>

    In order to import an STL file, use the following construct:
    \code
    STL { filename.stl }
    \endcode
    The \c .stl extension is not necessary, and the STL file can be in either
    ASCII or binary format. The file fomat specification can be found at
    http://en.wikipedia.org/wiki/STL_(file_format)

    The geometry must have consistent outward pointing normals.
    Problems will also occur in the geometry
    if there are any internal faces or vertices (that is, the mesh is not
		closed).
*/

/**
 * A class for the representation of boundaries in #D dimensions.
 *
 * \tparam D dimensionality of the space the boundary is embedded in
 * \tparam T numeric storage type
 */
template<int D, typename T>
class BoundaryMesh {
};

template<typename T>
class BoundaryMesh<2, T> {
public:
	typedef Vector<T, 2> Point;
	typedef Face<2, T> Face;
	typedef Vertex<2, T> Vertex;

private:
	std::vector<Face> faces_;
	std::vector<Vertex> vertices_;

public:
	BoundaryMesh<2, T>();

	BoundaryMesh<2, T>(const Polygon<T>& polygon);

	template<typename Grid>
	void scanConvert(Grid& grid) const;
};


template<typename T>
class BoundaryMesh<3, T> {
public:
	typedef Vector<T, 3> Point;
	typedef ConnectedFace<T> Face;
	typedef ConnectedEdge<T> Edge;
	typedef ConnectedVertex<T> Vertex;

	typedef std::vector<Face*> FaceContainer;
	typedef std::vector<Edge*> EdgeContainer;
	typedef std::vector<Vertex*> VertexContainer;

private:
	FaceContainer faces_;
	EdgeContainer edges_;
	VertexContainer vertices_;

public:
	//! Construct an empty boundary
	BoundaryMesh<3, T>() :
		faces_(0),
 		edges_(0),
		vertices_(0) { }

	/**
	 * Creates a BoundaryMesh from an STL file.
	 *
	 * Note that the STL format technically requires that the surface be closed.
	 * This is not fully tested for, but is assumed.  The surface need not be
	 * simply connected.
	 *
	 * The STL format also does not store adjacency relations; as such, vertices
	 * with the same coordinates are merged, and this information is used to
	 * generate face--face adjacency relations.  If coincident vertices are not
	 * exactly equal, this will fail.
	 *
	 * \param filename Path to the STL file.
	 */
	static BoundaryMesh<3, T> loadSTL(const std::string& filename);

	//@{
	// Accessors for vertices
	size_t vertices() const { return vertices_.size(); }
	typename VertexContainer::iterator verticesBegin() { return vertices_.begin(); }
	typename VertexContainer::iterator verticesEnd() { return vertices_.end(); }
	typename VertexContainer::const_iterator verticesBegin() const { return vertices_.begin(); }
	typename VertexContainer::const_iterator verticesEnd() const { return vertices_.end(); }
	//@}

	//@{
	// Accessors for edges
	size_t edges() const { return edges_.size(); }
	typename EdgeContainer::iterator edgesBegin() { return edges_.begin(); }
	typename EdgeContainer::iterator edgesEnd() { return edges_.end(); }
	typename EdgeContainer::const_iterator edgesBegin() const { return edges_.begin(); }
	typename EdgeContainer::const_iterator edgesEnd() const { return edges_.end(); }
	//@}

	//@{
	// Accessors for faces
	size_t faces() const { return faces_.size(); }
	typename FaceContainer::iterator facesBegin() { return faces_.begin(); }
	typename FaceContainer::iterator facesEnd() { return faces_.end(); }
	typename FaceContainer::const_iterator facesBegin() const { return faces_.begin(); }
	typename FaceContainer::const_iterator facesEnd() const { return faces_.end(); }
	//@}


	/** Fill grid with the signed distance function from the boundary.
	 *
	 * This uses the characteristic polyhedron/scan conversion algorithm described in
	 * Mauch, Sean "Efficient Algorithms for Solving Static Hamilton-Jacobi Equations" (2003) PhD Thesis, Caltech.
	 *
	 * This algorithm is linear in the number of points and the number of faces
	 * in the mesh O(N + F) and is thus optimal.
	 *
	 * \param grid The grid to fill with the signed distance field.
	 */
	template<typename Grid>
	void scanConvert(Grid& grid) const;
};

