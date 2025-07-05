#ifndef PROBLEMINFO
#define PROBLEMINFO

#include <Eigen/Dense>

struct ProblemInfo {
    double mu;
    double rho;
    double relax = 0.8;
    double relaxp = 0.8;
};

struct GridInfo {
    int NX;
    int NY;
    double LX;
    double LY;
    double dx;
    double dy;
    Eigen::VectorXd x;
    Eigen::VectorXd y;
};

// if a velocity boundary condition is prescribed: the pressure boundary condition at that wall is zero-gradient
// otherwise, a pressure boundary condition is used and the velocity boundary condition is left at zero-gradient
struct BoundaryConditions {
    // Top, Left, Right, Bottom
    // bool: flow field BC (or not)
    bool FIELD_T = true;
    bool FIELD_B = true;
    bool FIELD_L = true;
    bool FIELD_R = true;
    // U boundary conditions
    double U_T = 0;
    double U_L = 0;
    double U_R = 0;
    double U_B = 0;
    // V boundary conditions
    double V_T = 0;
    double V_L = 0;
    double V_R = 0;
    double V_B = 0;
    // P boundary conditions
    double P_T = 0;
    double P_L = 0;
    double P_R = 0;
    double P_B = 0;
};

#endif