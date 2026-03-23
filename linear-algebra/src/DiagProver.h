#ifndef DIAGPROVER_H
#define DIAGPROVER_H

#include <lapacke.h>
#include <cblas.h>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

class DiagProver
{
public:
    explicit DiagProver(const std::string &filepath);

    bool prove();

private:
    int n;
    std::vector<double> A;    // original matrix (row-major)
    std::vector<double> P;    // eigenvector matrix
    std::vector<double> D;    // diagonal eigenvalue matrix
    std::vector<double> PInv; // P inverse
    std::vector<double> eigvals, eigvals_im;

    // helpers
    bool isSymmetric() const;
    void eigenDecomposition(); // fills P, eigvals, eigvals_im
    void buildD();
    void buildPInverse();
    double frobeniusError() const;
    void printMatrix(const std::string &label, const std::vector<double> &M) const;
    void printEigenValues() const;

    static std::vector<double> matrixMultiplication(const std::vector<double> &A,
                                      const std::vector<double> &B, int n);
    static std::vector<double> luInvert(std::vector<double> M, int n);
};

#endif