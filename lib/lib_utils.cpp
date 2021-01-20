#include <vector>
#include "lib_utils.h"


/**
 * \brief Scales the data X into the unit hypercube
 *
 * sometimes data has a funny scale and it is outside of the libigl-viewer-frustrum
 * this is a convenience function to scale the data into the unit hypercube
 *
 * \params[in/out] X The data points organized as a n x d matrix
 */
void scale_to_unit_cube(Eigen::MatrixXd &X, double &s)
{
    if (s == 0)
    {
        s = 1.0 / X.colwise().maxCoeff().maxCoeff();
    }

    X = (X.array() * s).matrix();
}

/**
 * \brief Finds the length of the smallest edge in the mesh
 *
 * why? the convergence properties of the Laplacian depend on the smallest edge length
 *  so in case you need that you can scale your mesh accordingly
 *
 * \params[in] points The data points organized as a n x d matrix
 * \params[in] faces The faces (connectivity) organized as a n x d matrix
 */
double find_min_edge_length(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces)
{
    double min = std::numeric_limits<double>::max();

    for(int face_idx=0; face_idx < faces.rows(); face_idx++)
    {
        int i = faces(face_idx, 0);
        int j = faces(face_idx, 1);
        int k = faces(face_idx, 2);

        Eigen::Vector3d v_i = points.row(i);
        Eigen::Vector3d v_j = points.row(j);
        Eigen::Vector3d v_k = points.row(k);

        Eigen::Vector3d e_ij = v_j - v_i;
        Eigen::Vector3d e_jk = v_k - v_j;
        Eigen::Vector3d e_ki = v_i - v_k;

        min = std::min(e_ij.norm(), min);
        min = std::min(e_jk.norm(), min);
        min = std::min(e_ki.norm(), min);
    }

    return min;
}


void make_2d_grid(Eigen::Vector2d top_left, Eigen::Vector2d bottom_right, int size, Eigen::MatrixXd &grid_points, Eigen::MatrixXi &faces)
{
    Eigen::VectorXd x_lin_space = Eigen::VectorXd::LinSpaced(size, top_left(0), bottom_right(0));
    Eigen::VectorXd y_lin_space = Eigen::VectorXd::LinSpaced(size, top_left(1), bottom_right(1));

    grid_points.resize(size * size, 3);
    for(int i=0; i < size; i++)
    {
        for(int j=0; j < size; j++)
        {
            grid_points(i * size + j, 0) = x_lin_space(i);
            grid_points(i * size + j, 1) = y_lin_space(j);
            grid_points(i * size + j, 2) = 0.0;
        }
    }

    std::vector<Eigen::VectorXi> tris;
    for(int i=0; i < size - 1; i++) {
        for (int j = 0; j < size - 1; j++) {
            int curr = i * size + j;

            Eigen::Vector3i T1(curr, curr + 1, curr + 1 + size);
            Eigen::Vector3i T2(curr, curr + 1 + size, curr + size);

            tris.push_back(T1);
            tris.push_back(T2);
        }
    }

    faces.resize(tris.size(), 3);
    for(int i=0; i < tris.size(); i++)
    {
        faces.row(i) = tris[i];
    }
}