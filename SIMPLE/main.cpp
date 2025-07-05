#include <math.h>

#include <Eigen/Dense>
#include <iostream>

#include "InputReader.hpp"
#include "Boundary.hpp"
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

    // Boundary Conditions
    double U_LID = 2;

    // Fluid Properties
    double mu = 1;    // 0.0010518  # dynamic viscosity, Pa*s
    double rho = 10;  // 1000  # density, kg/m^3

    // Convergence
    double alpha = 0.5;
    double alpha_p = 0.5;

    ProblemInfo Problem;
    // Problem.mu = mu;
    // Problem.rho = rho;
    // Problem.relax = alpha;
    // Problem.relaxp = alpha_p;

    GridInfo Mesh;
    // Mesh.NX = NX;
    // Mesh.NY = NY;
    // Mesh.LX = LX;
    // Mesh.LY = LY;
    // Mesh.dx = dx;
    // Mesh.dy = dy;
    // Mesh.x = x;
    // Mesh.y = y;

    BoundaryConditions BC;
    // LID DRIVEN CAVITY
    //BC.U_T = U_LID;

    // DIFFERENTIAL PRESSURE PIPE FLOW?
    // TODO: not sure if the linear interp (extrapolation) used in
    //      the pressure boundary conditions cases to extend the interior of velocity fields is entirely correct
    // small perturbation to velocity fields is needed for the pressure driven problem to work
    // BC.FIELD_L = false;
    // BC.FIELD_R = false;
    // BC.P_L = 800;

    ReadInputFile(Problem, Mesh, BC);
    ProcessInputFile(Problem, Mesh, BC);

    cout << Problem.relax << endl;
    cout << Problem.relaxp << endl;
    cout << Problem.rho << endl;
    cout << Problem.mu << endl;
    cout << Mesh.NX << endl;
    cout << Mesh.NY << endl;
    cout << Mesh.LX << endl;
    cout << Mesh.LY << endl;
    cout << Mesh.dx << endl;
    cout << Mesh.dy << endl;

    SIMPLE(BC, Mesh, Problem);

    return 0;
}