#ifndef ATCG1_LIB_LEAST_SQUARES_CONFORMAL_MAPS_H
#define ATCG1_LIB_LEAST_SQUARES_CONFORMAL_MAPS_H

#include <Eigen/Core>

Eigen::MatrixXd least_squares_conformal_maps(const Eigen::MatrixXi &faces, const Eigen::MatrixXd &points);

#endif //ATCG1_LIB_LEAST_SQUARES_CONFORMAL_MAPS_H
