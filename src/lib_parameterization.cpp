#include "lib_parameterization.h"
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/find.h>
#include <igl/slice_into.h>

/**
 * \brief detect the boundary vertices for a given mesh
 *
 * TODO:
 *  -compute all the indices at the boundary of a mesh
 *  -make them unique
 *  -compute the boundary indices from the edges
 *  -compute the interior indices with igl::setdiff, you can use igl::colon to generate a list of all the indices
 *
 *  -plot the results for the different smoothing parameters
 *
 *  \param[in] points The points of a mesh (n x 3)
 *  \param[in] faces The faces of a mesh (m x 3)
 *  \param[in/out] boundary_indices The indices of the boundary (as vector)
 *  \param[in/out] interior_indices The indices of the interior (as vector)
 *  \param[in/out] boundary_edges The path of the boundary indices (as matrix)
 */
void detect_boundary_vertices(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, Eigen::VectorXi &boundary_indices, Eigen::VectorXi &interior_indices, Eigen::MatrixXi &boundary_edges)
{
}

/**
 * \brief calculate the length of the boundary path
 *
 * TODO:
 *  -iterate over boundary edges in the path and compute the lengths of the edges
 *
 *  -plot the results for the different smoothing parameters
 *
 *  \param[in] path The indices of the path along the boundary in consecutive order (m x 2)
 *  \param[in] points The points of a mesh (n x 3)
 *  \param[out] edge_lengths The lengths of the edges of the boundary (m x 1)
 */
Eigen::VectorXd path_length(const Eigen::MatrixXi &path, const Eigen::MatrixXd &points)
{
    Eigen::VectorXd edge_lengths(path.rows());
    return edge_lengths;
}

/**
 * \brief sort the boundary edges to a path
 *
 * TODO:
 *  -iterate over boundary edges and sort them into a matrix with first column = from index, second column = to index
 *
 *  \param[in] boundary_edges The boundary edges
 *  \param[out] path The path in consecutive order
 */
Eigen::MatrixXi detect_boundary_path(const Eigen::MatrixXi &boundary_edges)
{
    Eigen::MatrixXi boundary_path(boundary_edges.rows(), 2);
    return boundary_path;
}

/**
 * \brief map a path to the unit circle
 *
 * TODO:
 *  -each edge length in the path is a step around the unit circle
 *  -edge_lendths must sum to one (arc length)
 *  -select arbitrary start vertex
 *
 *  \param[in] edge_lengths The edge lengths
 *  \param[out] circle The 3D positions of the boundary vertices along a circle
 */
Eigen::MatrixXd map_boundary_edges_to_circle(const Eigen::VectorXd &edge_lengths)
{
    Eigen::MatrixXd circle(edge_lengths.rows(), 3);
    return circle;
}


//--------------------------------------------------------------------------

/**
 * \brief main convenience function for parameterization
 *
 * what it does:
 *  -calculate the weights for the selected method
 *  -construct the linear operator from the weights and respective indices
 *  -construct right hand side containing boundary conditions and constraints (fixed boundary)
 *  -solve for the interior unknown vertex positions in the unit disk
 *
 *  \param[in] faces The triangles of a mesh
 *  \param[in] points The points (vertices) of a mesh
 *  \param[in] method The method how to calculate the weights (uniform spring, chordal spring, wachspress, discrete harmonic, mean value)
 *  \param[in] boundary_condition The boundary condition (e.g. Dirichlet=0)
 *  \param[in] boundary_indices The boundary indices
 *  \param[in] boundary_constraints The fixed values for boundary vertices
 *  \param[out] uv The 3D positions in the unit circle ( one coordinate is 0)
 */
Eigen::MatrixXd parameterize_with_hard_boundary_constraints(const Eigen::MatrixXi &faces, const Eigen::MatrixXd &points, const Method &method, const Eigen::MatrixXd &boundary_condition, const Eigen::VectorXi &boundary_indices, const Eigen::MatrixXd &boundary_constraints)
{
    int nv = points.rows();
    std::vector<Eigen::Triplet<double>> coefficients = method_switcher(points, faces, method);

    Eigen::SparseMatrix<double> op = construct_operator(coefficients, boundary_indices, nv);

    Eigen::MatrixXd rhs = construct_rhs(boundary_condition, boundary_indices, boundary_constraints);

    Eigen::MatrixXd uv = solve(op, rhs);
    return uv;
}

/**
 * \brief construct the operator with hard constraints (replaced rows with deltas)
 *
 * TODO:
 *  -construct operator from coefficients (i, j, val)
 *  -construct the linear operator from the weights and respective indices
 *  -compute the rowsums: op * 1 = rowsums, where 1 is a vector of ones
 *  -substract the rowsums from the diagonal
 *  -the operator must be normalilzed
 *  -set the hard constraints by replacing rows with deltas
 *
 *  \param[in] coefficients The triplet vector containg the indices and the weights
 *  \param[in] boundary_indices The boundary indices
 *  \param[in] nv The number of vertices (for dimensions of operator)
 *  \param[out] op The sparse linear operator
 */
Eigen::SparseMatrix<double> construct_operator(const std::vector<Eigen::Triplet<double>> &coefficients, const Eigen::VectorXi &boundary_indices, int nv)
{
    Eigen::SparseMatrix<double> op(nv, nv);
    return op;
}

/**
 * \brief construct the right hand side of the system
 *
 * TODO:
 *  -replace the boundary condition at the boundary with the constraint
 *
 *  \param[in] boundary_condition The boundary condition (nv x 3)
 *  \param[in] boundary_indices The row indices (k x 1)
 *  \param[in] boundary_constraints The column indices (k x 3)
 *  \param[out] op The sparse linear operator
 */
Eigen::MatrixXd construct_rhs(const Eigen::MatrixXd &boundary_condition,
                              const Eigen::VectorXi &boundary_indices,
                              const Eigen::MatrixXd &boundary_constraints)
{
    Eigen::MatrixXd rhs = boundary_condition;
    return rhs;
}


/**
 * \brief Computes a numerically stable version of law of cosines for the angle from the metric
 *
 *  -You should use this function!
 *
 *  \param[in] a,b,c The edge lengths in unsorted order
 *  \param[out] alpha The angle between a and b, opposite to c
 */
inline double angle_from_metric(double a, double b, double c)
{
    /* numerically stable version of law of cosines
     * angle between a and b, opposite to edge c
    */

    double alpha = acos((a*a + b*b - c*c) / (2.0 * a * b));

    if (alpha < 1e-8)
    {
        alpha = std::sqrt((c*c - (a - b)*(a - b)) / (2.0 * a * b));
        std::cout << "small angle < 1e-8!" << std::endl;
    }
    return alpha;
}

void clip_and_warn(Eigen::MatrixXd &values)
{
    Eigen::VectorXd mean = values.cwiseAbs().colwise().mean();
    for(int i = 0; i < values.rows() - 1; i++)
    {
        for(int j=0; j < values.cols(); j++)
        {
            if (values(i, j) > mean(j))
            {
                values(i, j) = mean(j);
            }
            else if(values(i, j) < -mean(j))
            {
                values(i, j) = -mean(j);
            }
        }
    }
}

/**
 * \brief compute the weights and indices for the specified method
 *
 * TODO:
 *  -select the type of weights to fill the operator with
 *  -uniform
 *  -chordal
 *  -wachspress
 *  -discrete harmonic
 *  -mean value
 *
 *  HINT:
 *      - In a triangle with vertices i, j, k think of what values you have to set to wij and the other weights
 *      - Pay attention to the direction, is ij always equal to ji? (handle the other weights consistently)
 *      - Pick a consistent naming convention for vertices, edges, faces angles and lengths!
 *      - This will reduce a lot of errors!
 *
 *  \param[in] points The vertices of a mesh
 *  \param[in] faces The faces of a mesh
 *  \param[in] method The specified method for the computation of the weights
 *  \param[out] coefficients The filled triplet vector from which the operator is constructed
 */
std::vector<Eigen::Triplet<double>> method_switcher(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, const Method &method)
{
    int nt = faces.rows();
    std::string message;

    std::vector<Eigen::Triplet<double>> coefficients;

    for(int count = 0; count < nt; count ++)
    {
        int i = faces(count, 0);
        int j = faces(count, 1);
        int k = faces(count, 2);

        Eigen::Vector3d vi = points.row(i);
        Eigen::Vector3d vj = points.row(j);
        Eigen::Vector3d vk = points.row(k);

        Eigen::Vector3d eij = vj - vi;
        Eigen::Vector3d ejk = vk - vj;
        Eigen::Vector3d eki = vi - vk;

        double rij = eij.norm();
        double rjk = ejk.norm();
        double rki = eki.norm();

        double alphai = angle_from_metric(rij, rki, rjk);
        double alphaj = angle_from_metric(rjk, rij, rki);
        double alphak = angle_from_metric(rjk, rki, rij);

        double wij = 0, wjk = 0, wki = 0;
        double wji = 0, wkj = 0, wik = 0;
        //HINT: in a triangle with vertices i, j, k think of what values you have to set to wij and the other weights
        // Pay attention to the direction, is ij always equal to ji? (handle the other weights consistently)

        switch(method){
            case Method::UNIFORM_SPRING :
            {
                // implement here uniform weights
                message = "uniform spring";
            }break;

            case Method::CHORDAL_SPRING :
            {
                //implement here chordal spring weitghts: w = 1.0 / r^2
                message = "chordal spring";
            }break;

            case Method::WACHSPRESS :
            {
                //implement here the wachspress weights
                message = "wachspress";
            }break;

            case Method::DISCRETE_HARMONIC :
            {
                //implement here the discrete harmonic weights
                message = "discrete harmonic";
            }break;

            case Method::MEAN_VALUE :
            {
                //implement here the mean value weights
                message = "mean value";
            }break;
        };

        coefficients.emplace_back(i, j, wij);
        coefficients.emplace_back(j, k, wjk);
        coefficients.emplace_back(k, i, wki);

        //symmetric part
        coefficients.emplace_back(j, i, wji);
        coefficients.emplace_back(k, j, wkj);
        coefficients.emplace_back(i, k, wik);
    }


    std::cout << "compute parameterization with " << message << " weights" << std::endl;
    return coefficients;
}

/**
 * \brief solve the linear system
 *
 * what it does:
 *  -select which function to execute for the specified method
 *
 *  \param[in] op The operator (left hand side)
 *  \param[in] rhs The right hand side
 *  \param[out] x The solution of the linear system
 */
Eigen::MatrixXd solve(const Eigen::SparseMatrix<double> &op, const Eigen::MatrixXd &rhs)
{
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
    solver.compute(op);
    return solver.solve(rhs).eval();
}

/**
 * \brief compute the distortion measure for the parameterization
 *
 *  TODO:
 *      -implement the following distortion measures:
 *      -area distortion measure: s1 * s2 = 1
 *      -length distortion measure: s1 = s2 = 1
 *      -angle distortion measure: s1 / s2 = 1
 *      -bijectivity distortion measure: sign(s1 * s2)
 *      -conformalfactor distortion measure: s1
 *
 *  \param[in] parameterization The parameterization (flat)
 *  \param[in] points The vertices of a mesh
 *  \param[in] faces The triangles of a mesh
 *  \param[in] measure The specifier for a distortion measure
 *  \param[in/out] distortion The distortion value for each triangle in a vector
 */
void measure_distortion(const Eigen::MatrixXd &parameterization, const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, const Measure &measure, Eigen::VectorXd &distortion)
{
    Eigen::MatrixXd param(parameterization.rows(), 2);
    param.col(0) = parameterization.col(0);
    param.col(1) = parameterization.col(1);

    Eigen::MatrixXd q0s, q1s, q2s, qe0s, qe1s, v0s, v1s, v2s, ve0s, ve1s;
    Eigen::VectorXi all_points, all_0, all_1, all_2, all_01, all_012, all_faces, i0s, i1s, i2s;

    igl::colon<int>(0, faces.rows() - 1, all_faces);
    igl::colon<int>(0, param.rows() - 1, all_points);
    igl::colon<int>(0, 0, all_0);
    igl::colon<int>(1, 1, all_1);
    igl::colon<int>(2, 2, all_2);
    igl::colon<int>(0, 1, all_01);
    igl::colon<int>(0, 2, all_012);

    igl::slice(faces, all_faces, all_0, i0s);
    igl::slice(faces, all_faces, all_1, i1s);
    igl::slice(faces, all_faces, all_2, i2s);

    igl::slice(param, i0s, all_01, q0s);
    igl::slice(param, i1s, all_01, q1s);
    igl::slice(param, i2s, all_01, q2s);

    qe0s = q1s - q0s;
    qe1s = q2s - q0s;

    igl::slice(points, i0s, all_012, v0s);
    igl::slice(points, i1s, all_012, v1s);
    igl::slice(points, i2s, all_012, v2s);

    ve0s = v1s -v0s;
    ve1s = v2s -v0s;

    distortion.resize(faces.rows());

    for(int i=0; i < faces.rows(); i++)
    {
        Eigen::Matrix2d q(2, 2);
        q << qe0s.row(i), qe1s.row(i);

        Eigen::MatrixXd v(3, 2);
        v << ve0s.row(i).transpose(), ve1s.row(i).transpose();

        Eigen::MatrixXd A = v * q.inverse();
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

        //S are the singular values from 2d to 3d
        Eigen::VectorXd S = svd.singularValues();

        switch(measure){
            case Measure::AREA :
            {
                //implement here: area disortion measure
            }break;

            case Measure::ANGLE :
            {
                //implement here: angle distortion measure
            }break;

            case Measure::LENGTH:
            {
                //implemnt here: length distortion measure
            }break;

            case Measure::BIJECTIVITY :
            {
                //implement here: bijectivity measure (is there a flip or not?)
            }break;

            case Measure::CONFORMALFACTOR :
            {
                //implement here: conformal factor
            }break;
        };
    }
    std::cout << "distortion max: " << distortion.maxCoeff() << " min: " << distortion.minCoeff() << std::endl;
}


