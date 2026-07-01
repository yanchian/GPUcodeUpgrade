
//! Initialise from a polygon, assuming that the points are oriented clockwise
template<typename T>
BoundaryMesh<2, T>::
BoundaryMesh<2, T>(const Polygon<T>& polygon) {
	for (int i = 0; i < polygon.points(); i++) {
		Point left   = polygon.point(i - 1),
		      middle = polygon.point(i),
		      right  = polygon.point(i + 1);
		Point leftDisplacement = middle - left;
		Point rightDisplacement = right - middle;
		Point leftNormal(-leftDisplacement[1], leftDisplacement[0]);
		Point rightNormal(-rightDisplacement[1], rightDisplacement[0]);
		leftNormal /= abs(leftNormal);
		rightNormal /= abs(rightNormal);

		// initialise a face with the left points
		faces_.push_back(Face(left, middle, leftNormal));

		// initialise a vertex using all the data
		vertices_.push_back(Vertex(middle, leftNormal, rightNormal));
	}
}

template<typename T>
template<typename Grid>
void
BoundaryMesh<2, T>::
scanConvert(Grid& grid) const {
//	const Grid2DWriter<Grid, T> writer(grid, distance_);
	const T range = 2.0 * std::max(grid.xMax() - grid.xMin(), grid.yMax() - grid.yMin());
	for (typename std::vector<Face>::const_iterator iter = faces_.begin(); iter != faces_.end(); ++iter) {
		(*iter).createScannableCharacteristicPolygon(range).scanConvert(grid);
		(*iter).createScannableCharacteristicPolygon(-range).scanConvert(grid);
	}
	for (typename std::vector<Vertex>::const_iterator iter = vertices_.begin(); iter != vertices_.end(); ++iter) {
		(*iter).createScannableCharacteristicPolygon(range).scanConvert(grid);
		(*iter).createScannableCharacteristicPolygon(-range).scanConvert(grid);
	}
}

template<typename T>
struct VertexTemp {
	ConnectedVertex<T>* vertex;
	size_t id;
};

template<typename T>
bool operator<(const VertexTemp<T>& a, const VertexTemp<T>& b) {
	return *(a.vertex) < *(b.vertex);
}

template<typename T>
struct EdgeTemp {
	ConnectedEdge<T>* edge;
	size_t id;
};

template<typename T>
bool operator<(const EdgeTemp<T>& a, const EdgeTemp<T>& b) {
	return *(a.edge) < *(b.edge);
}

template<typename T>
BoundaryMesh<3, T>
BoundaryMesh<3, T>::
loadSTL(const std::string& filename) {
	/*#if GLOBAL_DIMENSION != 3
	  Log::log("STL geometry only available in three dimensions", Log::UserAbort);
	#else*/
	BoundaryMesh<3, T> mesh;
	std::ifstream stlFile;

	stlFile.open(filename.c_str(), std::ios_base::in);

	if (!stlFile) {
		// throw Error::Conversion() << Error::Text("Unable to open file " + filename);
	}

	char header[80];

	// Take first 5 bytes: solid => ASCII o'wise assume binary
	stlFile.get(header, 6);

	std::vector<VertexTemp<T> > vertexTemp;

	if (strncmp(header, "solid", 5)  == 0 || strncmp(header, "SOLID", 5) == 0) {
		// ASCII file
		// Get solid name - ignored
		std::string solidname;
		std::getline(stlFile, solidname);

		try {
			while (!stlFile.eof()) {
				std::string s1;
				std::string s2;
				Face* f = new Face();
				// facet normal
				stlFile >> s1 >> s2;
				// n1 n2 n3
				stlFile >> f->normal()[0] >> f->normal()[1] >> f->normal()[2];

				// If the file has ended, s1, s2 = "endsolid" "name of solid"
				// and then getting n1,n2,n3 will have failed
				if (!stlFile) {
					delete f;
					break;
				}

				// outer loop
				stlFile >> s1 >> s2;
				// vertex v11 v12 v13
				for (int i = 0 ; i < 3 ; i++) {
					Vertex* v = new Vertex();
					stlFile >> s1;
					for (int j = 0 ; j < 3 ; j++) {
						stlFile >> v->position()[j];
					}
					VertexTemp<T> vt;
				 	vt.vertex = v;
					vt.id = vertexTemp.size();
					vertexTemp.push_back(vt);

					f->vertex(i) = (Vertex*) vt.id;
				}
				// endloop
				stlFile >> s1;
				// endfacet
				stlFile >> s1;
				mesh.faces_.push_back(f);
			}
		} catch (std::ios_base::failure) {
			//throw Error::Conversion() << Error::Text("Unexpected EOF");
			throw;
		}

	} else {
		stlFile.close();
		// reopen in binary format
		stlFile.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
		stlFile.read(header, 80);

		// Take care that we use an unsigned 32-bit integer
		boost::uint32_t facets;
		stlFile.read((char*)&facets, 4);
		// STL spec. requires 4-byte floating point number
		// with IEEE 754 format
		if (sizeof(float) != 4) {
			//Log::log("TriangulatedSurface requires a 4 byte floating-point number, and float isn't suitable", Log::DeveloperAbort);
		}

		try {
			for (boost::uint32_t i = 0; i < facets ; i++) {
				Face* f = new Face();
				Point normal;
				for (int j = 0; j < 3; j++) {
					float tmp;
					stlFile.read((char*)&tmp, 4);
					normal[j] = tmp;
				}
				f->normal(normal);

				for (int j = 0 ; j < 3 ; j++) {
					Vertex* v = new Vertex();
					Point point;
					for (int k = 0 ; k < 3 ; k++) {
						float tmp;
						stlFile.read((char*)&tmp, 4);
						point[k] = tmp;
					}
					v->position(point);

					VertexTemp<T> vt;
				 	vt.vertex = v;
					vt.id = vertexTemp.size();
					vertexTemp.push_back(vt);

					// I am going to hell for this
					f->vertex(j) = (Vertex*) vt.id;
				}
				mesh.faces_.push_back(f);

				boost::uint16_t byteCount;
				stlFile.read((char*)&byteCount, 2); // unused
			}
		} catch (std::ios_base::failure) {
			//throw Error::Conversion() << Error::Text("Unexpected EOF");
			throw;
		}
	}
	stlFile.close();

	// sort the vertices so that identical ones are adjacent
	std::sort(vertexTemp.begin(), vertexTemp.end());

	// nuke the duplicates
	// construct a vector which maps id -> Vertex* with duplicates all removed
	std::vector<Vertex*> idVectorMap(vertexTemp.size());
	VertexTemp<T> lastVT;
	lastVT.vertex = NULL;
	for (typename std::vector<VertexTemp<T> >::iterator iter = vertexTemp.begin(); iter != vertexTemp.end(); ++iter) {
		if (lastVT.vertex != NULL && *(lastVT.vertex) == *(iter->vertex)) {
			delete iter->vertex;
			iter->vertex = lastVT.vertex;
		} else {
			mesh.vertices_.push_back(iter->vertex);
		}
		idVectorMap[iter->id] = iter->vertex;
		lastVT = *iter;
	}


	// go through all the faces and fix up their vertices with real pointers
	// also construct a vector of edges (with duplicates):
	std::vector<EdgeTemp<T> > edgeTemp;
	for (typename FaceContainer::iterator fIter = mesh.faces_.begin(); fIter != mesh.faces_.end(); ++fIter) {
		for (typename Face::VertexContainer::iterator vIter = (*fIter)->verticesBegin(); vIter != (*fIter)->verticesEnd(); ++vIter) {
			*vIter = idVectorMap[(size_t)(*vIter)];
			(*vIter)->facesPush(*fIter);
		}
		// append all the edges of this triangle
		for (int i = 0; i < 3; i++) {
			Edge* e = new Edge;
			e->vertex(0) = (*fIter)->vertex(i);
			e->vertex(1) = (*fIter)->vertex((i + 1) % 3);
			EdgeTemp<T> et = {e, edgeTemp.size()};
			edgeTemp.push_back(et); // should optimize this: we know the size in advance
			// going to hell for this, too
			(*fIter)->edge(i) = (Edge*) et.id;
		}
	}

	std::sort(edgeTemp.begin(), edgeTemp.end());
	// nuke the duplicates
	// construct a vector which maps id -> Edge* with duplicates all removed
	std::vector<Edge*> idEdgeMap(edgeTemp.size());
	EdgeTemp<T> lastET;
	lastET.edge = NULL;
	for (typename std::vector<EdgeTemp<T> >::iterator iter = edgeTemp.begin(); iter != edgeTemp.end(); ++iter) {
		if (lastET.edge != NULL && *(lastET.edge) == *(iter->edge)) {
			delete iter->edge;
			iter->edge = lastET.edge;
		} else {
			mesh.edges_.push_back(iter->edge);
		}
		idEdgeMap[iter->id] = iter->edge;
		lastET = *iter;
	}

	// go through all the faces and fix up their vertices with real pointers
	for (typename FaceContainer::iterator fIter = mesh.faces_.begin(); fIter != mesh.faces_.end(); ++fIter) {
		for (typename Face::EdgeContainer::iterator eIter = (*fIter)->edgesBegin(); eIter != (*fIter)->edgesEnd(); ++eIter) {
			*eIter = idEdgeMap[(size_t)*eIter];
			(*eIter)->facesPushBack(*fIter);
		}
		(*fIter)->init();
	}

	// assign edges to vertices
	for (typename EdgeContainer::iterator eIter = mesh.edges_.begin(); eIter != mesh.edges_.end(); ++eIter) {
		for (int i = 0; i < 2; i++) {
			(*eIter)->vertex(i)->edgePush(*eIter);
		}
		(*eIter)->init();
	}

	for (typename VertexContainer::iterator vIter = mesh.vertices_.begin(); vIter != mesh.vertices_.end(); ++vIter) {
		(*vIter)->init();
	}

/*#endif*/
	return mesh;
}

template<typename T>
template<typename Grid>
void
BoundaryMesh<3, T>::
scanConvert(Grid& grid) const {
	const T range = 2.0 * std::max(std::max(grid.xMax() - grid.xMin(), grid.yMax() - grid.yMin()), grid.zMax() - grid.zMin());
	for (typename FaceContainer::const_iterator iter = faces_.begin(); iter != faces_.end(); ++iter) {
		std::cout << **iter << std::endl;
		(*iter)->createScannableCharacteristicPolyhedron(range).scanConvert(grid);
		(*iter)->createScannableCharacteristicPolyhedron(-range).scanConvert(grid);
	}
	for (typename EdgeContainer::const_iterator iter = edges_.begin(); iter != edges_.end(); ++iter) {
		std::cout << **iter << std::endl;
		if ((*iter)->convex()) {
			(*iter)->createScannableCharacteristicPolyhedron(range).scanConvert(grid);
		}
		if ((*iter)->concave()) {
			(*iter)->createScannableCharacteristicPolyhedron(-range).scanConvert(grid);
		}
	}
	for (typename VertexContainer::const_iterator iter = vertices_.begin(); iter != vertices_.end(); ++iter) {
		std::cout << **iter << std::endl;
		if ((*iter)->convex()) {
				(*iter)->createScannableCharacteristicPolyhedron(range).scanConvert(grid);
		}
		if ((*iter)->concave()) {
			(*iter)->createScannableCharacteristicPolyhedron(-range).scanConvert(grid);
		}
	}
}

template<typename T>
std::ostream&
operator<<(std::ostream& out, const BoundaryMesh<3, T>& mesh) {
	out << "BoundaryMesh<3> {" << std::endl;
	out << "\tfaces (" << mesh.faces() << " ): [" << std::endl;
	for (typename BoundaryMesh<3, T>::FaceContainer::const_iterator iter =  mesh.facesBegin(); iter != mesh.facesEnd(); ++iter) {
		out << "\t\t" << **iter << "," << std::endl;
	}
	out << "\t]," << std::endl << "\tvertices(" << mesh.vertices() << "): [" << std::endl;
	for (typename BoundaryMesh<3, T>::VertexContainer::const_iterator iter =  mesh.verticesBegin(); iter != mesh.verticesEnd(); ++iter) {
		out << "\t\t" << **iter << ", " << std::endl;
	}
	out << "\t]," << std::endl << "\tedges(" << mesh.edges() << "): [" << std::endl;
	for (typename BoundaryMesh<3, T>::EdgeContainer::const_iterator iter =  mesh.edgesBegin(); iter != mesh.edgesEnd(); ++iter) {
		out << "\t\t" << **iter << ", " << std::endl;
	}
	out << "\t]" << std::endl;
	out << "}";
	return out;
}


