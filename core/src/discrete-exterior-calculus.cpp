// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    const size_t number_of_vertices = mesh.getVertexIndices().size();

    auto hodge_star_matrix = Eigen::SparseMatrix<double>(number_of_vertices, number_of_vertices);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(number_of_vertices);
    for(size_t i = 0; i < number_of_vertices; i++){
        const auto vertex = mesh.vertex(i); 
        //for each face the dual vertex is situated in such a way, that it's quadrilateral's 
        //area is exactly 1/3 of face's area
        const double area_of_dual_face = barycentricDualArea(vertex);
        triplets.push_back({i, i, area_of_dual_face});
    }
    hodge_star_matrix.setFromTriplets(triplets.begin(), triplets.end());

    return hodge_star_matrix;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    const size_t number_of_edges = mesh.getEdgeIndices().size();

    SparseMatrix<double> hodge_star_matrix(number_of_edges, number_of_edges);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(number_of_edges);
    for(size_t i = 0; i < number_of_edges; i++){
        const auto& edge = mesh.edge(i);
        const auto& half_edge = edge.halfedge();
        const auto volume_relation = 0.5 * (cotan(half_edge) + cotan(half_edge.twin()));
        triplets.push_back({i, i, volume_relation});
    }
    hodge_star_matrix.setFromTriplets(triplets.begin(), triplets.end());
    return hodge_star_matrix;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    const size_t number_of_faces = mesh.getFaceIndices().size();

    auto hodge_star_matrix = SparseMatrix<double>(number_of_faces, number_of_faces);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(number_of_faces);
    for(size_t i = 0; i < number_of_faces; i++){
        const auto face = mesh.face(i);
        const auto volume_relation = 1.0 / faceArea(face);
        triplets.push_back({i, i, volume_relation});
    }
    hodge_star_matrix.setFromTriplets(triplets.begin(), triplets.end());

    return hodge_star_matrix;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
    const size_t number_of_vertices = mesh.getVertexIndices().size();
    const size_t number_of_edges = mesh.getEdgeIndices().size();

    auto derivative_matrix = SparseMatrix<double>(number_of_edges, number_of_vertices);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(number_of_edges * 2);
    for(size_t i = 0; i < number_of_edges; i++){
        const auto edge = mesh.edge(i);
        const auto halfedge = edge.halfedge();
        size_t tip_index = halfedge.tipVertex().getIndex();
        size_t tail_index = halfedge.tailVertex().getIndex();
        triplets.push_back({i, tip_index, 1.0});
        triplets.push_back({i, tail_index, -1.0});
    }

    derivative_matrix.setFromTriplets(triplets.begin(), triplets.end());

    return derivative_matrix;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    const size_t number_of_edges = mesh.getEdgeIndices().size();
    const size_t number_of_faces = mesh.getFaceIndices().size();

    auto derivative_matrix = SparseMatrix<double>(number_of_faces, number_of_edges);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(number_of_faces * 3);
    for(size_t i = 0; i < number_of_faces; i++){
        const auto face = mesh.face(i);
        const auto counter_clock_half_edges = face.adjacentHalfedges();
        for(const auto edge : face.adjacentEdges()){
            const size_t edge_index = edge.getIndex();
            const auto halfedge = edge.halfedge();
            const auto it = std::find(counter_clock_half_edges.begin(), counter_clock_half_edges.end(), halfedge);
            const double relative_rotation = it == counter_clock_half_edges.end() ? -1.0 : 1.0;
            triplets.push_back({i, edge_index, relative_rotation});
        }
    }

    derivative_matrix.setFromTriplets(triplets.begin(), triplets.end());

    return derivative_matrix;
}

} // namespace surface
} // namespace geometrycentral
