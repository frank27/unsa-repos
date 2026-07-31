#include "pca_lapack.h"
#include <cmath>
#include <algorithm>
#include <lapacke.h>

PCA::PCA(int num_components) : num_components(num_components) {}

std::vector<std::vector<double>> PCA::fit_transform(const std::vector<std::vector<double>>& data) {
    if (data.empty() || data[0].empty()) return {};

    int n_samples = data.size();
    int n_features = data[0].size();
    int k = std::min(num_components, n_features);

    // 1. Center the data (subtract mean from each feature)
    std::vector<double> mean(n_features, 0.0);
    for (const auto& row : data) {
        for (int j = 0; j < n_features; ++j) {
            mean[j] += row[j];
        }
    }
    for (int j = 0; j < n_features; ++j) {
        mean[j] /= n_samples;
    }

    std::vector<std::vector<double>> centered(n_samples, std::vector<double>(n_features));
    for (int i = 0; i < n_samples; ++i) {
        for (int j = 0; j < n_features; ++j) {
            centered[i][j] = data[i][j] - mean[j];
        }
    }

    // 2. Compute covariance matrix (n_features x n_features)
    // C = (X^T * X) / (n_samples - 1)
    std::vector<double> cov(n_features * n_features, 0.0);
    double scale = 1.0 / (n_samples - 1);
    for (int i = 0; i < n_features; ++i) {
        for (int j = 0; j <= i; ++j) {
            double sum = 0.0;
            for (int s = 0; s < n_samples; ++s) {
                sum += centered[s][i] * centered[s][j];
            }
            cov[i * n_features + j] = sum * scale;
            cov[j * n_features + i] = cov[i * n_features + j];  // Symmetric
        }
    }

    // 3. Compute eigenvalues and eigenvectors using LAPACKE dsyev
    std::vector<double> eigenvalues(n_features);

    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'L', n_features,
                              cov.data(), n_features, eigenvalues.data());

    if (info != 0) {
        return {};  // Error in LAPACK
    }

    // 4. Sort eigenvalues in descending order and select top k
    // LAPACK returns eigenvalues in ascending order, so we reverse
    std::vector<int> indices(n_features);
    for (int i = 0; i < n_features; ++i) indices[i] = i;
    std::sort(indices.begin(), indices.end(),
              [&eigenvalues](int a, int b) { return eigenvalues[a] > eigenvalues[b]; });

    // 5. Project data onto top k principal components
    // cov now contains eigenvectors in columns (after dsyev with "V")
    std::vector<std::vector<double>> result(n_samples, std::vector<double>(k));
    for (int s = 0; s < n_samples; ++s) {
        for (int comp = 0; comp < k; ++comp) {
            double proj = 0.0;
            int eig_idx = indices[comp];
            for (int f = 0; f < n_features; ++f) {
                proj += centered[s][f] * cov[f * n_features + eig_idx];
            }
            result[s][comp] = proj;
        }
    }

    return result;
}
