#include <cstdlib>
#include <iostream>

//set this to supress libigl viewer help
#define IGL_VIEWER_VIEWER_QUIET

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/parula.h>
#include <igl/jet.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/colon.h>

#include <Eigen/Core>
#include "lib_utils.h"
#include "lib_parameterization.h"
#include "lib_least_squares_conformal_maps.h"

/**
 * \brief First exercise
 *
 * what it does:
 *  -visualize a convex polygon
 *  -visualize the weights by varying the position of the interior point on a surface
 *
 * TODO:
 *  -complete the method_switcher function by implementing the missing weights
 *  -test each weight type (uniform, chordal, wachspress, discrete harmonic, mean value)
 *      by changing the parameter in switcher method: Method::UNIFORM_SPRING etc..
 *
 *  -plot the results for the different parameterization parameters
 *
 */
void exercise1()
{
    Eigen::MatrixXd V(7, 3);
    V <<    2, 2, 0,
            6, 0, 0,
            6, 4, 0,
            2.5, 4, 0,
            0, 3, 0,
            1, 2, 0,
            0, 1, 0;

    Eigen::MatrixXi F(6, 3);
    F << 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 5, 0, 5, 6, 0, 6, 1; // not convex

    //Eigen::MatrixXi F(5, 3);
    //F << 0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 6, 0, 6, 1; // convex

    int idxMiddle = 0;
    V = V.rowwise() - V.colwise().mean();
    Eigen::RowVector3d middle = V.row(idxMiddle);

    Eigen::VectorXi boundary_indices, interior_indices;
    Eigen::MatrixXi boundary_edges;
    detect_boundary_vertices(V, F, boundary_indices, interior_indices, boundary_edges);

    int size = 100;

    double _min = V.minCoeff();
    double _max = V.maxCoeff();

    Eigen::MatrixXd grid_points;
    Eigen::MatrixXi grid_faces;

    make_2d_grid(Eigen::Vector2d(_min - 1, _min - 1), Eigen::Vector2d(_max + 1, _max + 1), size, grid_points, grid_faces);

    int nv = V.rows();
    Eigen::MatrixXd scalars(grid_points.rows(), nv);

    for(int i=0; i < grid_points.rows(); i++)
    {
        V.row(idxMiddle) = grid_points.row(i);
        std::vector<Eigen::Triplet<double>> coefficients = method_switcher(V, F, Method::DISCRETE_HARMONIC);

        Eigen::SparseMatrix<double> op(nv, nv);
        op.setFromTriplets(coefficients.begin(), coefficients.end());

        scalars.row(i) = op.row(0);
    }

    clip_and_warn(scalars);

    Eigen::MatrixXi T(F.rows() * 3, 2);
    T << F.col(0), F.col(1), F.col(1), F.col(2), F.col(2), F.col(0);

    Eigen::MatrixXd colors(scalars.rows(), 3);


    V.row(idxMiddle) = middle;

    for(int i=0; i < nv; i++)
    {
        Eigen::VectorXd s = scalars.col(i);
        if(s.cwiseAbs().mean() < 1e-9)
        {
            continue;
        }

        igl::jet(s, true, colors);

        igl::opengl::glfw::Viewer viewer;
        viewer.data().set_mesh(grid_points, grid_faces);
        viewer.data().set_colors(colors);

        viewer.data().set_edges(V, T, Eigen::RowVector3d(0.0, 0.0, 0.0));
        viewer.launch();
    }
}

/**
 * \brief Second exercise
 *
 * what it does:
 *  -load a mesh with boundary
 *  -detect the path along the boundary of the mesh
 *  -separate vertex indices into interior and boundary indices
 *  -map the boundary to the unit disk
 *  -compute the parameterization for uniform weights
 *  -compute the parameterization for chordal weights
 *  -compute the parameterization for wachspress weights
 *  -compute the parameterization for discrete harmonic weights
 *  -compute the parameterization for mean value weights
 *  -compute the parameterization for conformal weights (least squares conformal maps, lscm)
 *  -select a distortion measure (area, length, angle, bijectivity, conformality)
 *
 *  -plot the results for the different parameterization parameters
 *
 *  \param[in] filename1 The path to the mesh
 */
void exercise2(const std::string &filename1){
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    bool load_success_target = igl::readOFF(filename1, V, F);
    if (!load_success_target)
    {
        std::cerr << "could not load file: " << filename1 << std::endl;
        return;
    }

    double s = 0.0;
    scale_to_unit_cube(V, s);
    V = V.rowwise() - V.colwise().mean();

    Measure measure = Measure::ANGLE;

    Eigen::VectorXi boundary_indices, interior_indices, all;
    Eigen::MatrixXi boundary_edges;
    Eigen::MatrixXd boundary_condition = Eigen::MatrixXd::Zero(V.rows(), 3);

    detect_boundary_vertices(V, F, boundary_indices, interior_indices, boundary_edges);
    Eigen::MatrixXi path = detect_boundary_path(boundary_edges);

    boundary_indices = path.col(0);

    Eigen::VectorXd edge_lengths = path_length(path, V);
    Eigen::MatrixXd circle = map_boundary_edges_to_circle(edge_lengths);

    Eigen::VectorXd distortion;
    Eigen::MatrixXd colors(F.rows(), 3);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.launch();

    Eigen::MatrixXd uv_uniform_spring = parameterize_with_hard_boundary_constraints(F, V, Method::UNIFORM_SPRING, boundary_condition, boundary_indices, circle);

    measure_distortion(uv_uniform_spring, V, F, measure, distortion);
    igl::jet(distortion, true, colors);

    igl::opengl::glfw::Viewer viewer1;
    viewer1.data().set_mesh(uv_uniform_spring, F);
    viewer1.data().set_colors(colors);
    viewer1.data().add_points(circle, Eigen::RowVector3d(0.0, 1.0, 0.0));
    viewer1.data().point_size = 3.0;
    viewer1.launch();

    Eigen::MatrixXd uv_chordal_spring = parameterize_with_hard_boundary_constraints(F, V, Method::CHORDAL_SPRING, boundary_condition, boundary_indices, circle);

    measure_distortion(uv_chordal_spring, V, F, measure, distortion);
    igl::jet(distortion, true, colors);

    igl::opengl::glfw::Viewer viewer2;
    viewer2.data().set_mesh(uv_chordal_spring, F);
    viewer2.data().set_colors(colors);
    viewer2.data().add_points(circle, Eigen::RowVector3d(0.0, 1.0, 0.0));
    viewer2.data().point_size = 3.0;
    viewer2.launch();

    Eigen::MatrixXd uv_wachspress = parameterize_with_hard_boundary_constraints(F, V, Method::WACHSPRESS, boundary_condition, boundary_indices, circle);

    measure_distortion(uv_wachspress, V, F, measure, distortion);
    igl::jet(distortion, true, colors);

    igl::opengl::glfw::Viewer viewer3;
    viewer3.data().set_mesh(uv_wachspress, F);
    viewer3.data().set_colors(colors);
    viewer3.data().add_points(circle, Eigen::RowVector3d(0.0, 1.0, 0.0));
    viewer3.data().point_size = 3.0;
    viewer3.launch();

    Eigen::MatrixXd uv_discrete_harmonic = parameterize_with_hard_boundary_constraints(F, V, Method::DISCRETE_HARMONIC, boundary_condition, boundary_indices, circle);

    measure_distortion(uv_discrete_harmonic, V, F, measure, distortion);
    igl::jet(distortion, true, colors);

    igl::opengl::glfw::Viewer viewer4;
    viewer4.data().set_mesh(uv_discrete_harmonic, F);
    viewer4.data().set_colors(colors);
    viewer4.data().add_points(circle, Eigen::RowVector3d(0.0, 1.0, 0.0));
    viewer4.data().point_size = 3.0;
    viewer4.launch();

    Eigen::MatrixXd uv_mean_value = parameterize_with_hard_boundary_constraints(F, V, Method::MEAN_VALUE, boundary_condition, boundary_indices, circle);

    measure_distortion(uv_mean_value, V, F, measure, distortion);
    igl::jet(distortion, true, colors);

    igl::opengl::glfw::Viewer viewer5;
    viewer5.data().set_mesh(uv_mean_value, F);
    viewer5.data().set_colors(colors);
    viewer5.data().add_points(circle, Eigen::RowVector3d(0.0, 1.0, 0.0));
    viewer5.data().point_size = 3.0;
    viewer5.launch();

    Eigen::MatrixXd uv_lscm = least_squares_conformal_maps(F, V);

    measure_distortion(uv_lscm, V, F, measure, distortion);
    igl::jet(distortion, true, colors);

    igl::opengl::glfw::Viewer viewer6;
    viewer6.data().set_mesh(uv_lscm, F);
    viewer6.data().set_colors(colors);
    viewer6.launch();
    }


/**
* \brief The main function called when running this program
*
* what it does:
*  -check provided filenames
*  -run both exercises in a row
*
*  \param[in] argc The number of arguments to the binary
*  \param[in] argv The array of arguments to the binary
*/
int main(int argc, char *argv[])
{
    std::string filename1, filename2;

    if (argc == 2)
    {
        filename1 = argv[1];
    }
    else
    {
        std::cerr << "please call assignmentsheet7 like this: " << std::endl;
        std::cerr << "./bin/assignmentsheet7 data/maxear.off" << std::endl;
        return EXIT_FAILURE;
    }
    //exercise1();
    exercise2(filename1);
	return EXIT_SUCCESS;
}
