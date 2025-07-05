#include <Eigen/Dense>
#include <iostream>
#include <sstream>
#include <string>

#include "ProblemInfo.hpp"

using namespace std;

void ReadInputFile(ProblemInfo& Problem, GridInfo& Mesh, BoundaryConditions& BC) {
    string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) {
            continue;
        }

        stringstream ss(line);
        string string_field;
        double double_field;

        if (ss >> string_field >> double_field) {
            if (string_field == string("mu")) {
                Problem.mu = double_field;
            } else if (string_field == string("rho")) {
                Problem.rho = double_field;
            } else if (string_field == string("relax")) {
                Problem.relax = double_field;
            } else if (string_field == string("relaxp")) {
                Problem.relaxp = double_field;
            } else if (string_field == string("NX")) {
                Mesh.NX = (int)double_field;
            } else if (string_field == string("NY")) {
                Mesh.NY = (int)double_field;
            } else if (string_field == string("LX")) {
                Mesh.LX = double_field;
            } else if (string_field == string("LY")) {
                Mesh.LY = double_field;
            } else if (string_field == string("FIELD_T")) {
                BC.FIELD_T = (bool)double_field;
            } else if (string_field == string("FIELD_L")) {
                BC.FIELD_L = (bool)double_field;
            } else if (string_field == string("FIELD_R")) {
                BC.FIELD_R = (bool)double_field;
            } else if (string_field == string("FIELD_B")) {
                BC.FIELD_B = (bool)double_field;
            } else if (string_field == string("U_T")) {
                BC.U_T = double_field;
            } else if (string_field == string("U_L")) {
                BC.U_L = double_field;
            } else if (string_field == string("U_R")) {
                BC.U_R = double_field;
            } else if (string_field == string("U_B")) {
                BC.U_B = double_field;
            } else if (string_field == string("V_T")) {
                BC.V_T = double_field;
            } else if (string_field == string("V_L")) {
                BC.V_L = double_field;
            } else if (string_field == string("V_R")) {
                BC.V_R = double_field;
            } else if (string_field == string("V_B")) {
                BC.V_B = double_field;
            } else if (string_field == string("P_T")) {
                BC.P_T = double_field;
            } else if (string_field == string("P_L")) {
                BC.P_L = double_field;
            } else if (string_field == string("P_R")) {
                BC.P_R = double_field;
            } else if (string_field == string("P_B")) {
                BC.P_B = double_field;
            } else {
            }
        } else {
            cerr << "Warning: Could not parse line: \""
                 << line << "\"" << endl;
        }
    }
}

void ProcessInputFile(ProblemInfo& Problem, GridInfo& Mesh, BoundaryConditions& BC) {
    Mesh.dx = Mesh.LX / (Mesh.NX - 1);
    Mesh.dy = Mesh.LY / (Mesh.NY - 1);

    // x and y locations
    Mesh.x = Eigen::VectorXd::LinSpaced(Mesh.NX, 0, Mesh.LX);
    Mesh.y = Eigen::VectorXd::LinSpaced(Mesh.NY, 0, Mesh.LY);
}