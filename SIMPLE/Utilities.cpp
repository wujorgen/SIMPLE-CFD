#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using Eigen::ArrayXXd;
using namespace std;

void WriteArrayXXdToCSV(const ArrayXXd& arr, const string& fname) {
    ofstream file(fname);

    for (int i = 0; i < arr.rows(); i++) {
        for (int j = 0; j < arr.cols(); j++) {
            file << arr(i, j);
            if (j < arr.cols() - 1)
                file << ", ";
        }
        file << endl;
    }

    file.close();
}