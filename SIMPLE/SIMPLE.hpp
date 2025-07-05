#ifndef SIMPLE_HPP
#define SIMPLE_HPP

#ifndef PROBLEMINFO
#include "ProblemInfo.hpp"
#endif

void SIMPLE(const BoundaryConditions&, const GridInfo&, const ProblemInfo&);

void FinalUVPFields(
    Eigen::ArrayXXd&,
    Eigen::ArrayXXd&,
    Eigen::ArrayXXd&,
    const Eigen::ArrayXXd&,
    const Eigen::ArrayXXd&,
    const Eigen::ArrayXXd&);

void EigenRef(Eigen::ArrayXXd&);

#endif