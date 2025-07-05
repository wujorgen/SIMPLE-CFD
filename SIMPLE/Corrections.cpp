#include <Eigen/Dense>

#include "ProblemInfo.hpp"

using Eigen::ArrayXXd;

void CalcPressureCorrection(ArrayXXd& p_corr, ArrayXXd& p_b, const ArrayXXd& u_star, const ArrayXXd& v_star, const ArrayXXd& d_e, const ArrayXXd& d_n, const GridInfo& Mesh) {
    double a_E;
    double a_W;
    double a_N;
    double a_S;
    double a_P;

    for (int i = 1; i < p_corr.rows() - 1; i++) {
        for (int j = 1; j < p_corr.cols() - 1; j++) {
            p_b(i, j) = -(u_star(i, j) - u_star(i, j - 1)) * Mesh.dy - (v_star(i - 1, j) - v_star(i, j)) * Mesh.dx;
        }
    }
    for (int g = 0; g < 10; g++) { // TODO
        for (int i = 1; i < p_corr.rows() - 1; i++) {
            for (int j = 1; j < p_corr.cols() - 1; j++) {
                a_E = -d_e(i, j) * Mesh.dy;
                a_W = -d_e(i, j - 1) * Mesh.dy;
                a_N = -d_n(i - 1, j) * Mesh.dx;
                a_S = -d_n(i, j) * Mesh.dx;

                a_P = a_E + a_W + a_N + a_S;

                p_corr(i, j) = a_E * p_corr(i, j + 1) + a_W * p_corr(i, j - 1) + a_N * p_corr(i - 1, j) + a_S * p_corr(i + 1, j);
                p_corr(i, j) += p_b(i, j);
                p_corr(i, j) /= a_P;
            }
        }
    }
}

void ApplyPressureCorrection(ArrayXXd& p, const ArrayXXd& p_corr, const ProblemInfo& Problem) {
    for (int i = 1; i < p.rows() - 1; i++) {
        for (int j = 1; j < p.cols() - 1; j++) {
            p(i, j) += Problem.relaxp * p_corr(i, j);
        }
    }
}

void ApplyUCorrection(ArrayXXd& u, const ArrayXXd& u_star, const ArrayXXd& d_e, const ArrayXXd& p_corr, const ProblemInfo& Problem) {
    for (int i = 1; i < u.rows() - 1; i++) {
        for (int j = 1; j < u.cols() - 1; j++) {
            u(i, j) = u_star(i, j) + Problem.relax * d_e(i, j) * (p_corr(i, j + 1) - p_corr(i, j));
        }
    }
}

void ApplyVCorrection(ArrayXXd& v, const ArrayXXd& v_star, const ArrayXXd& d_n, const ArrayXXd& p_corr, const ProblemInfo& Problem) {
    for (int i = 1; i < v.rows() - 1; i++) {
        for (int j = 1; j < v.cols() - 1; j++) {
            v(i, j) = v_star(i, j) + Problem.relax * d_n(i, j) * (p_corr(i, j) - p_corr(i + 1, j));
        }
    }
}