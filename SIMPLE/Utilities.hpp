#ifndef UTILITIES
#define UTILITIES

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template <typename T>
void WriteVectorToCSV(const std::vector<T>& vec, const std::string& fname) {
    std::ofstream file(fname);

    // file << "itr, error" << std::endl;

    for (size_t i = 0; i < vec.size(); i++) {
        file << vec[i] << std::endl;
        // file << i << ", " << vec[i] << std::endl;
    }

    file.close();
};

void WriteArrayXXdToCSV(const Eigen::ArrayXXd&, const std::string&);

#endif