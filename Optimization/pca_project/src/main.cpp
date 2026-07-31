#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "pca_lapack.h"

std::vector<std::vector<double>> read_csv(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            try {
                row.push_back(std::stod(cell));
            } catch (...) {
                continue;
            }
        }
        if (!row.empty()) data.push_back(row);
    }
    return data;
}

int main() {
    std::cout << "=== PCA with LAPACK Example ===" << std::endl;

    auto data = read_csv(DATA_DIR "wine.csv");
    if (data.empty()) {
        std::cerr << "Failed to load data" << std::endl;
        return 1;
    }

    std::cout << "Loaded " << data.size() << " samples with " << data[0].size() << " features" << std::endl;

    int n_components = 2;
    PCA pca(n_components);

    auto transformed = pca.fit_transform(data);

    std::cout << "\nTransformed to " << n_components << " components:" << std::endl;
    std::cout << "First 5 samples:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), transformed.size()); ++i) {
        std::cout << "Sample " << i + 1 << ": [";
        for (size_t j = 0; j < transformed[i].size(); ++j) {
            std::cout << std::fixed << std::setprecision(4) << transformed[i][j];
            if (j < transformed[i].size() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }

    return 0;
}
