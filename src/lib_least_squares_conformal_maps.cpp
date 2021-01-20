#include <iostream>
#include "lib_least_squares_conformal_maps.h"
#include <Eigen/Sparse>
#include <igl/colon.h>
#include <igl/slice.h>

/**
 * \brief least squares conformal maps
 *
 * TODO:
 *  -iterate over the faces
 *  -compute the edges
 *  -compute local frame from edges
 *  -project edges into local frame
 *  -construct the gradient maps
 *  -collect the values and indices in a coefficients vector
 *
 *  -build the linear system
 *  -fix one edge (two vertices)
 *  -solve for the unknown vertex positions in the plane
 *
 *  \param[in] faces The triangles of a mesh
 *  \param[in] points The vertices of a mesh
 *  \param[out] uv The parameterization (coordinates in the plane (3D with one coordinate being 0))
 */
Eigen::MatrixXd least_squares_conformal_maps(const Eigen::MatrixXi &faces,
                                             const Eigen::MatrixXd &points)
{
    std::cout << "compute parameterization with least squares conformal maps (lscm)" << std::endl;
    Eigen::MatrixXd uv(points.rows(), 2);

    return uv;
}
