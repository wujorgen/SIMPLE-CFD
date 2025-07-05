#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "Boundary.hpp"
#include "Corrections.hpp"
#include "MomentumEquations.hpp"
#include "ProblemInfo.hpp"
#include "Utilities.hpp"

using Eigen::all;
using Eigen::ArrayXXd;
using Eigen::last;
using Eigen::seq;
using Eigen::VectorXd;
using namespace std;

void EigenRef(ArrayXXd& testarr) {
    // order matters. even tho this function is declared in header,
    //      it must be defined before use.
    testarr(0, all) = 1234;
}

void FinalUVPFields(ArrayXXd& u_final, ArrayXXd& v_final, ArrayXXd& p_final,
                    const ArrayXXd& u, const ArrayXXd& v, const ArrayXXd& p) {
    for (int i = 0; i < p.rows() - 1; i++) {
        for (int j = 0; j < p.cols() - 1; j++) {
            // p grid is offset in both directions
            p_final(i, j) = (p(i, j) + p(i + 1, j) + p(i, j + 1) + p(i + 1, j + 1)) / 4;
            // u grid is offset in y direction / along i index
            u_final(i, j) = (u(i, j) + u(i + 1, j)) / 2;
            // v grid is offsetin x direction / along j index
            v_final(i, j) = (v(i, j) + v(i, j + 1)) / 2;
        }
    }
}

void SIMPLE(const BoundaryConditions& BC, const GridInfo& Mesh, const ProblemInfo& Problem) {
    // PRESSURE ghost grids: (NY + 1, NX + 1)
    // The pressure ghost grid is fully staggered from the vertices.
    ArrayXXd p = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX + 1);
    ArrayXXd p_star = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX + 1);
    ArrayXXd p_corr = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX + 1);
    ArrayXXd p_b = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX + 1);
    // p.fill(0.001);
    ApplyPBoundary(p, BC);

    // u grid: (NY + 1, NX)
    // staggered only in y-direction
    // note that row 0 (top) and row NY + 1 (bottom) must be zero, as those are ghost velocities.
    ArrayXXd u = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX);
    ArrayXXd u_star = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX);
    ArrayXXd u_corr = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX);
    ArrayXXd d_e = ArrayXXd::Zero(Mesh.NY + 1, Mesh.NX);
    u.fill(0.001);
    ApplyUBoundary(u, BC);

    // v grid: (NY, NX + 1)
    // staggered only in x-direction
    // note that col 0 (left) and col NX + 1 (right) must be zero, as those are ghost velocities.
    ArrayXXd v = ArrayXXd::Zero(Mesh.NY, Mesh.NX + 1);
    ArrayXXd v_star = ArrayXXd::Zero(Mesh.NY, Mesh.NX + 1);
    ArrayXXd v_corr = ArrayXXd::Zero(Mesh.NY, Mesh.NX + 1);
    ArrayXXd d_n = ArrayXXd::Zero(Mesh.NY, Mesh.NX + 1);
    v.fill(0.001);
    ApplyVBoundary(v, BC);

    int itr = 0;
    int minitr = 10;
    int maxitr = 1000;
    double error = 1;
    double u_error;
    double v_error;
    double ethresh = 1e-5;
    bool DIVERGED = false;
    bool CONVERGED = false;

    vector<double> errors;

    // TODO: timestep loop around this, with conditional to disable transient term
    while (!DIVERGED && (error > ethresh || itr < minitr) && (itr < maxitr)) {
        // calc u-momentum
        CalcUStar(u_star, u, v, p, d_e, Mesh, Problem);

        // apply u-momentum boundary conditions
        ApplyUBoundary(u_star, BC);

        // call v-momentum
        CalcVStar(v_star, u, v, p, d_n, Mesh, Problem);

        // apply v-momentum boundary conditions
        ApplyVBoundary(v_star, BC);

        // zero out pressure corrections
        p_corr.setZero();
        p_b.setZero();

        // calculate pressure corrections
        CalcPressureCorrection(p_corr, p_b, u_star, v_star, d_e, d_n, Mesh);  // cout << p_b << endl;

        // apply pressure corrections to pressure field
        ApplyPressureCorrection(p, p_corr, Problem);

        // apply pressure boundary conditions
        // NoPressureGradientAtBoundary(p);
        ApplyPBoundary(p, BC);

        // correct u-momentum
        ApplyUCorrection(u, u_star, d_e, p_corr, Problem);

        // apply u-momentum boundary conditions
        ApplyUBoundary(u, BC);

        // correct v-momentum
        ApplyVCorrection(v, v_star, d_n, p_corr, Problem);

        // apply v-momentum boundary conditions
        ApplyVBoundary(v, BC);

        // check for convergence
        error = p_b.matrix().norm();
        // error = max((u - u_star).maxCoeff(), (v - v_star).maxCoeff());
        errors.push_back(error);
        if (error > 1) {
            cout << (error > 1) << endl;
            cout << "SOLUTION DIVERGED" << endl;
            DIVERGED = true;
        }
        cout << "ITR: " << itr << ", ERROR: " << error << endl;
        itr++;
    }

    if (itr <= maxitr && error <= ethresh) {
        CONVERGED = true;
        cout << "SOLUTION CONVERGED" << endl;
    }

    ArrayXXd u_final = ArrayXXd::Zero(Mesh.NY, Mesh.NX);
    ArrayXXd v_final = ArrayXXd::Zero(Mesh.NY, Mesh.NX);
    ArrayXXd p_final = ArrayXXd::Zero(Mesh.NY, Mesh.NX);

    FinalUVPFields(u_final, v_final, p_final, u, v, p);

    WriteVectorToCSV(errors, string("convergencehistory.txt"));

    WriteArrayXXdToCSV(p_final, string("p.txt"));
    WriteArrayXXdToCSV(u_final, string("u.txt"));
    WriteArrayXXdToCSV(v_final, string("v.txt"));

    std::vector<double> x_vec(Mesh.x.data(), Mesh.x.data() + Mesh.x.size());
    std::vector<double> y_vec(Mesh.y.data(), Mesh.y.data() + Mesh.y.size());

    WriteVectorToCSV(x_vec, string("x_domain.txt"));
    WriteVectorToCSV(y_vec, string("y_domain.txt"));

    // cout << p << endl;
    // EigenRef(p);
    // cout << p << endl;
    // cout << u_star << endl;
}
