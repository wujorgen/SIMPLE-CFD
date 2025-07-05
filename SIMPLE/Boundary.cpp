#include <Eigen/Dense>
#include <iostream>
#include "ProblemInfo.hpp"

using Eigen::all;
using Eigen::ArrayXXd;
using Eigen::last;
using Eigen::seq;
using Eigen::VectorXd;
using namespace std;

void ApplyUBoundary(ArrayXXd &u, const BoundaryConditions &BC) {
    // left
    if (BC.FIELD_L) {
        u(seq(1, last - 1), 0) = BC.U_L;
    } else {
        u(seq(1, last - 1), 0) = 2 * u(seq(1, last - 1), 1) - u(seq(1, last - 1), 2);
    }
    // right
    if (BC.FIELD_R) {
        u(seq(1, last - 1), last) = BC.U_R;
    } else {
        u(seq(1, last - 1), last) = 2 * u(seq(1, last - 1), last - 1) - u(seq(1, last - 1), last - 2);
    }
    // top
    if (BC.FIELD_T) {
        u(0, all) = BC.U_T * 2 - u(1, all);
    } else {
        u(0, all) = u(1, all);
    }
    // bottom
    if (BC.FIELD_B) {
        u(last, all) = BC.U_B * 2 - u(last - 1, all);
    } else {
        u(last, all) = u(last - 1, all);
    }
}

void ApplyVBoundary(ArrayXXd &v, const BoundaryConditions &BC) {
    // top
    if (BC.FIELD_T) {
        v(0, seq(1, last - 1)) = BC.V_T;
    } else {
        v(0, seq(1, last - 1)) = 2 * v(1, seq(1, last - 1)) - v(2, seq(1, last - 1));
    }
    // bottom
    if (BC.FIELD_B) {
        v(last, seq(1, last - 1)) = BC.V_B;
    } else {
        v(last, seq(1, last - 1)) = 2 * v(last - 1, seq(1, last - 1)) - v(last - 2, seq(1, last - 1));
    }
    // left
    if (BC.FIELD_L) {
        v(all, 0) = BC.V_L * 2 - v(all, 1);
    } else {
        v(all, 0) = v(all, 1);
    }
    // right
    if (BC.FIELD_R) {
        v(all, last) = BC.V_R * 2 - v(all, last - 1);
    } else {
        v(all, last) = v(all, last - 1);
    }
}

void NoPressureGradientAtBoundary(ArrayXXd &p) {
    // top
    p(0, all) = p(1, all);
    // bottom
    p(last, all) = p(last - 1, all);
    // left
    p(all, 0) = p(all, 1);
    // right
    p(all, last) = p(all, last - 1);
}

void ApplyPBoundary(ArrayXXd &p, const BoundaryConditions &BC) {
    // top
    if (BC.FIELD_T) {
        p(0, all) = p(1, all);
    } else {
        //p(0, all) = BC.P_T * 2 - p(1, all);
        p(0, all) = BC.P_T;
        p(1, all) = BC.P_T;
    }
    // bottom
    if (BC.FIELD_B) {
        p(last, all) = p(last - 1, all);
    } else {
        // p(last, all) = BC.P_B * 2 - p(last - 1, all);
        p(last, all) = BC.P_B;
        p(last - 1, all) = BC.P_B;
    }
    // left
    if (BC.FIELD_L) {
        p(all, 0) = p(all, 1);
    } else {
        // p(all, 0) = BC.P_L * 2 - p(all, 1);
        p(all, 0) = BC.P_L;
        p(all, 1) = BC.P_L;
    }
    // right
    if (BC.FIELD_R) {
        p(all, last) = p(all, last - 1);
    } else {
        // p(all, last) = BC.P_R * 2 - p(all, last - 1);
        p(all, last) = BC.P_R;
        p(all, last - 1) = BC.P_R;
    }
}