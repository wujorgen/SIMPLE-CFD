#include <math.h>

#include <Eigen/Dense>
#include <iostream>

#include "Boundary.hpp"
#include "InputReader.hpp"
#include "SIMPLE.hpp"

using Eigen::all;
using Eigen::ArrayXXd;
using Eigen::last;
using Eigen::seq;
using Eigen::VectorXd;
using namespace std;

int main() {
    int NX = 51;
    double LX = 0.5;

    int NY = 31;
    double LY = 0.3;

    double dx = LX / (NX - 1);
    double dy = LY / (NY - 1);

    // x and y locations
    VectorXd x = VectorXd::LinSpaced(NX, 0, LX);
    VectorXd y = VectorXd::LinSpaced(NY, 0, LY);

    // note that the matrices are indexed as:
    // 0-------j
    // |
    // |
    // i

    // while the grid is indexed as
    // Y
    // |
    // |
    // 0-------X

    // The grid defines vertices.

    ProblemInfo Problem;

    GridInfo Mesh;

    BoundaryConditions BC;

    ReadInputFile(Problem, Mesh, BC);

    ProcessInputFile(Problem, Mesh, BC);

    SIMPLE(BC, Mesh, Problem);

    return 0;
}