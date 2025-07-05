#ifndef CORRECTIONS
#define CORRECTIONS

#include <Eigen/Dense>

#include "ProblemInfo.hpp"

void CalcPressureCorrection(
    Eigen::ArrayXXd&,        // pressure corrections
    Eigen::ArrayXXd&,        // pressure source (?)
    const Eigen::ArrayXXd&,  // u_star
    const Eigen::ArrayXXd&,  // v_star
    const Eigen::ArrayXXd&,  // d_e
    const Eigen::ArrayXXd&,  // d_n
    const GridInfo&          //
);

void ApplyPressureCorrection(
    Eigen::ArrayXXd&,        // pressure field
    const Eigen::ArrayXXd&,  // pressure corrections
    const ProblemInfo&);

void ApplyUCorrection(
    Eigen::ArrayXXd&,
    const Eigen::ArrayXXd&,
    const Eigen::ArrayXXd&,
    const Eigen::ArrayXXd&,
    const ProblemInfo&);

void ApplyVCorrection(
    Eigen::ArrayXXd&,
    const Eigen::ArrayXXd&,
    const Eigen::ArrayXXd&,
    const Eigen::ArrayXXd&,
    const ProblemInfo&);

#endif