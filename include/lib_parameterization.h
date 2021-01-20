#ifndef ATCG7_LIB_PARAMETERIZATION1_CPP_H
#define ATCG7_LIB_PARAMETERIZATION1_CPP_H
#include <Eigen/Core>
#include <Eigen/Sparse>

enum class Measure{
    AREA,
    LENGTH,
    ANGLE,
    BIJECTIVITY,
    CONFORMALFACTOR
};

enum class Method{
    UNIFORM_SPRING,
    CHORDAL_SPRING,
    WACHSPRESS,
    DISCRETE_HARMONIC,
    MEAN_VALUE
};

void detect_boundary_vertices(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, Eigen::VectorXi &boundary_indices, Eigen::VectorXi &interior_indices, Eigen::MatrixXi &boundary_edges);
Eigen::MatrixXi detect_boundary_path(const Eigen::MatrixXi &boundary_edges);
Eigen::VectorXd path_length(const Eigen::MatrixXi &path, const Eigen::MatrixXd &points);
Eigen::MatrixXd map_boundary_edges_to_circle(const Eigen::VectorXd &edge_lengths);

Eigen::MatrixXd parameterize_with_hard_boundary_constraints(const Eigen::MatrixXi &faces,
                                                            const Eigen::MatrixXd &points,
                                                            const Method &method,
                                                            const Eigen::MatrixXd &boundary_condition,
                                                            const Eigen::VectorXi &boundary_indices,
                                                            const Eigen::MatrixXd &boundary_constraints);

std::vector<Eigen::Triplet<double>> method_switcher(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, const Method &method);

Eigen::SparseMatrix<double> construct_operator(const std::vector<Eigen::Triplet<double>> &coefficients, const Eigen::VectorXi &boundary_indices, int nv);

Eigen::MatrixXd construct_rhs(const Eigen::MatrixXd &boundary_condition,
                              const Eigen::VectorXi &boundary_indices,
                              const Eigen::MatrixXd &boundary_constraints);

Eigen::MatrixXd solve(const Eigen::SparseMatrix<double> &op,
                      const Eigen::MatrixXd &rhs);

void measure_distortion(const Eigen::MatrixXd &parameterization,
                        const Eigen::MatrixXd &points,
                        const Eigen::MatrixXi &faces,
                        const Measure &measure,
                        Eigen::VectorXd &distortion);

void clip_and_warn(Eigen::MatrixXd &values);
#endif //ATCG7_LIB_PARAMETERIZATION1_CPP_H
