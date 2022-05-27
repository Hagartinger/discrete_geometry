// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
	using T = Eigen::Triplet<size_t>;

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

	std::vector<T> triplets;
	for(const auto& e : mesh->edges()){
		const auto edge_index = e.getIndex();
		// 1 stands for 'present'
		const auto value = 1;
		triplets.push_back({edge_index, e.firstVertex().getIndex(), value});
		triplets.push_back({edge_index, e.secondVertex().getIndex(), value});
	}

	auto result = SparseMatrix<size_t>(geometry->edgeIndices.size(), geometry->vertexIndices.size());
	result.setFromTriplets(triplets.begin(), triplets.end());

    return result;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
	using T = Eigen::Triplet<size_t>;

    geometry->requireEdgeIndices();
	geometry->requireFaceIndices();

	std::vector<T> triplets;

	for(const auto& face : mesh->faces()){
		for( const auto& edge : face.adjacentEdges()){
			triplets.push_back({face.getIndex(), edge.getIndex(), 1});
		}
	}

	auto result = SparseMatrix<size_t>(geometry->faceIndices.size(), geometry->edgeIndices.size());
	result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
	geometry->requireVertexIndices();

	auto result = Vector<size_t>(geometry->vertexIndices.size());
	for( auto& value : result){
		value = 0;
	}
	for( const auto& e : subset.vertices){
		result(e) = 1;
	}

    return result;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
	geometry->requireEdgeIndices();

	auto result = Vector<size_t>(geometry->edgeIndices.size());
	for( auto& value : result){
		value = 0;
	}
	for( const auto& e : subset.edges){
		result(e) = 1;
	}

    return result;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
	geometry->requireFaceIndices();

	auto result = Vector<size_t>(geometry->faceIndices.size());
	for( auto& value : result){
		value = 0;
	}
	for( const auto& e : subset.faces){
		result(e) = 1;
	}

    return result;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */

MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
	geometry->requireVertexIndices();
	geometry->requireEdgeIndices();
	geometry->requireFaceIndices();

	auto result = subset;
	for(const auto& sub_index: subset.vertices){
		const auto& vertex = mesh->vertex(sub_index);
		for(const auto& edge : vertex.adjacentEdges()){
			result.addEdge(edge.getIndex());
		}
		for(const auto& face : vertex.adjacentFaces()){
			result.addFace(face.getIndex());
		}
	}
	
	for(const auto& sub_index: subset.edges){
		const auto& edge = mesh->edge(sub_index);
		for( const auto& face : edge.adjacentFaces()){
			result.addFace(face.getIndex());
		}
	}
    return result;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
	auto result = subset;
	for(const auto& sub_index: subset.edges){
		const auto& edge = mesh->edge(sub_index);
		for(const auto& vertex : edge.adjacentVertices()){
			result.addVertex(vertex.getIndex());
		}
	}
	for(const auto& sub_index: subset.faces){
		const auto& face = mesh->face(sub_index);
		for(const auto& vertex : face.adjacentVertices()){
			result.addVertex(vertex.getIndex());
		}
		for(const auto& edge : face.adjacentEdges()){
			result.addEdge(edge.getIndex());
		}
	}
    return result;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
	auto result = closure(star(subset));
	result.deleteSubset(star(closure(subset)));
    return result;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */

bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
	for(const auto face_index : subset.faces){
		const auto face = mesh->face(face_index);
		for(const auto& vertex : face.adjacentVertices()){
			if(!subset.vertices.contains(vertex.getIndex())){
				return false;
			}
		}
		for(const auto& edge : face.adjacentEdges()){
			if(!subset.edges.contains(edge.getIndex())){
				return false;
			}
		}
	}

	for(const auto edge_index : subset.edges){
		const auto edge = mesh->edge(edge_index);
		for(const auto& vertex : edge.adjacentVertices()){
			if(!subset.vertices.contains(vertex.getIndex())){
				return false;
			}
		}
	}
    return true; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */

bool is_empty(const MeshSubset& set){
	return set.vertices.empty() && 
		   set.edges.empty() && 
		   set.faces.empty();
}

int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

	MeshSubset face_simplex;
	for(const auto& face_index : subset.faces){
		const auto face = mesh->face(face_index);
		face_simplex.addFace(face_index);
		for(const auto& edge : face.adjacentEdges()){
			face_simplex.addEdge(edge.getIndex());
		}
		for(const auto& vert : face.adjacentVertices()){
			face_simplex.addVertex(vert.getIndex());
		}
	}
	if(!is_empty(face_simplex)){
		if( face_simplex.equals(subset) ){
			return 2;
		} else {
			return -1;
		}
	}

	MeshSubset edge_simplex;
	for(const auto& edge_index : subset.edges){
		edge_simplex.addEdge(edge_index);
		const auto& edge = mesh->edge(edge_index);
		for(const auto& vertex : edge.adjacentVertices()){
			edge_simplex.addVertex(vertex.getIndex());
		}
	}

	if(!is_empty(edge_simplex)){
		if(edge_simplex.equals(subset)){
			return 1;
		} else {
			return -1;
		}
	}

	if(!is_empty(subset)){
		return 0;
	} else {
		return -1;
	}
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
	MeshSubset result;

	const auto degree = isPureComplex(subset);
	if( degree == 2){
		for(const auto edge_index : subset.edges){
			const auto& edge = mesh->edge(edge_index);
			bool is_boundary = false;
			for(const auto& face : edge.adjacentFaces()){
				if(!subset.faces.contains(face.getIndex())){
					is_boundary = true;
					break;
				}
			}
			if(is_boundary) {
				result.addEdge(edge_index);
				for(const auto& vert : edge.adjacentVertices()){
					result.addVertex(vert.getIndex());
				}
			}
		}
	} else if(degree == 1){
		for(const auto vert_index : subset.vertices){
			const auto& vertex = mesh->vertex(vert_index);
			size_t subset_degree = 0;
			for(const auto& e : vertex.adjacentEdges()){
				if(subset.edges.contains(e.getIndex())){
					subset_degree++;
				}
			}
			if(subset_degree == 1) {
				result.addVertex(vert_index);
			}
		}
	}
    return result; // placeholder
}
