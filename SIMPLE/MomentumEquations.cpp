#include <math.h>

#include <Eigen/Dense>
#include <iostream>

#include "ProblemInfo.hpp"

using Eigen::all;
using Eigen::ArrayXXd;
using Eigen::last;
using Eigen::seq;
using Eigen::VectorXd;
using namespace std;

void CalcUStar(ArrayXXd& u_star, ArrayXXd& u, ArrayXXd& v, ArrayXXd& p, ArrayXXd& d_e, const GridInfo Mesh, ProblemInfo Problem) {
    // velocities at faces of u-velocity control volume
    double u_E;
    double u_W;
    double v_N;
    double v_S;
    // u-momentum equation coefficients
    double a_E;
    double a_W;
    double a_N;
    double a_S;
    double a_e;
    double A_e;

    for (int i = 1; i < u_star.rows() - 1; i++) {
        for (int j = 1; j < u_star.cols() - 1; j++) {
            // interp @ faces of u-velocity control volume
            u_E = (u(i, j) + u(i, j + 1)) / 2;
            u_W = (u(i, j) + u(i, j - 1)) / 2;
            // v velocity grid is offset relative to u-velocity grid
            v_N = (v(i - 1, j) + v(i - 1, j + 1)) / 2;
            v_S = (v(i, j) + v(i, j + 1)) / 2;
            // u-momo eqn coeff
            a_E = -(u_E / 2) * Mesh.dy + (Problem.mu / Problem.rho) * (Mesh.dy / Mesh.dx);
            a_W = (u_W / 2) * Mesh.dy + (Problem.mu / Problem.rho) * (Mesh.dy / Mesh.dx);
            a_N = -(v_N / 2) * Mesh.dx + (Problem.mu / Problem.rho) * (Mesh.dx / Mesh.dy);
            a_S = (v_S / 2) * Mesh.dx + (Problem.mu / Problem.rho) * (Mesh.dx / Mesh.dy);

            a_e = (u_E / 2) * Mesh.dy - (u_W / 2) * Mesh.dy + (v_N / 2) * Mesh.dx - (v_S / 2) * Mesh.dx;
            a_e += (Problem.mu / Problem.rho) * ((1 / Mesh.dy + 1 / Mesh.dy) * Mesh.dx + (1 / Mesh.dx + 1 / Mesh.dx) * Mesh.dy);

            A_e = -Mesh.dy / Problem.rho;

            // store d_e for presure correction
            d_e(i, j) = A_e / a_e;

            u_star(i, j) = (a_E * u(i, j + 1) + a_W * u(i, j - 1) + a_N * u(i - 1, j) + a_S * u(i + 1, j)) / a_e + d_e(i, j) * (p(i, j + 1) - p(i, j));
        }
    }
}

void CalcVStar(ArrayXXd& v_star, ArrayXXd& u, ArrayXXd& v, ArrayXXd& p, ArrayXXd& d_n, const GridInfo Mesh, ProblemInfo Problem) {
    // velocities at faces of u-velocity control volume
    double u_E;
    double u_W;
    double v_N;
    double v_S;
    // u-momentum equation coefficients
    double a_E;
    double a_W;
    double a_N;
    double a_S;
    double a_n;
    double A_n;

    for (int i = 1; i < v_star.rows() - 1; i++) {
        for (int j = 1; j < v_star.cols() - 1; j++) {
            // interp @ faces of v-velocity control volume
            u_E = (u(i, j) + u(i + 1, j)) / 2;
            u_W = (u(i, j - 1) + u(i + 1, j - 1)) / 2;
            // u velocity grid is offset relative to v-velocity grid
            v_N = (v(i - 1, j) + v(i, j)) / 2;
            v_S = (v(i + 1, j) + v(i, j)) / 2;

            // v-momo eqn coeff
            a_E = -(u_E / 2) * Mesh.dy + (Problem.mu / Problem.rho) * (Mesh.dy / Mesh.dx);
            a_W = (u_W / 2) * Mesh.dy + (Problem.mu / Problem.rho) * (Mesh.dy / Mesh.dx);
            a_N = -(v_N / 2) * Mesh.dx + (Problem.mu / Problem.rho) * (Mesh.dx / Mesh.dy);
            a_S = (v_S / 2) * Mesh.dx + (Problem.mu / Problem.rho) * (Mesh.dx / Mesh.dy);

            a_n = (u_E / 2) * Mesh.dy - (u_W / 2) * Mesh.dy + (v_N / 2) * Mesh.dx - (v_S / 2) * Mesh.dx;
            a_n += (Problem.mu / Problem.rho) * ((1 / Mesh.dx + 1 / Mesh.dx) * Mesh.dy + (1 / Mesh.dy + 1 / Mesh.dy) * Mesh.dx);

            A_n = -Mesh.dx / Problem.rho;

            // store d_n for pressure correction
            d_n(i, j) = A_n / a_n;

            v_star(i, j) = (a_E * v(i, j + 1) + a_W * v(i, j - 1) + a_N * v(i - 1, j) + a_S * v(i + 1, j)) / a_n + d_n(i, j) * (p(i, j) - p(i + 1, j));
        }
    }
}