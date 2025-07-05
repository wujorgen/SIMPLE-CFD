#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include <Eigen/Dense>
#include "ProblemInfo.hpp"

void ApplyUBoundary(Eigen::ArrayXXd&, const BoundaryConditions&);

void ApplyVBoundary(Eigen::ArrayXXd&, const BoundaryConditions&);

void NoPressureGradientAtBoundary(Eigen::ArrayXXd&);

void ApplyPBoundary(Eigen::ArrayXXd&, const BoundaryConditions&);

#endif