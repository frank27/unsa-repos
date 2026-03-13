#ifndef STANDARDALGORITHM_H
#define STANDARDALGORITHM_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "Matrix.h"

class StandardAlgorithm
{
    Matrix multiplyStandard(const Matrix &A, const Matrix &B, int n)
    {
        Matrix C(n, std::vector<double>(n, 0.0));

        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k)
            {
                double aik = A[i][k];
                for (int j = 0; j < n; ++j)
                    C[i][j] += aik * B[k][j];
            }

        return C;
    }

public:
    Matrix multiply(const Matrix &A, const Matrix &B, int n)
    {
        return multiplyStandard(A, B, n);
        // Utils::print(a, n, "Standard Algorithm");
        // return a;
    }
};

#endif