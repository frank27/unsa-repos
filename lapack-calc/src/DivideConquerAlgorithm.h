#ifndef DIVIDECONQUERALGORITHM_H
#define DIVIDECONQUERALGORITHM_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "Matrix.h"

class DivideConquerAlgorithm
{
    static constexpr int CUTOFF = 64;

public:
    using Matrix = std::vector<std::vector<double>>;

    /* Public entry point */
    Matrix multiply(const Matrix &A, const Matrix &B, int n)
    {
        return multiplyDivideConquer(A, B, n);
        // Utils::print(a, n, "Divide & Conquer");
        // return a;
    }

private:
    /* Standard O(n^3) multiplication */
    Matrix multiplyStandard(const Matrix &A, const Matrix &B)
    {
        int n = A.size();
        Matrix C(n, std::vector<double>(n, 0.0));

        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++)
                for (int j = 0; j < n; j++)
                    C[i][j] += A[i][k] * B[k][j];

        return C;
    }

    /* Matrix addition */
    Matrix add(const Matrix &A, const Matrix &B)
    {
        int n = A.size();
        Matrix C(n, std::vector<double>(n));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                C[i][j] = A[i][j] + B[i][j];

        return C;
    }

    /* Divide & Conquer recursive algorithm */
    Matrix multiplyDivideConquer(const Matrix &A, const Matrix &B, int n)
    {
        /* Base case */
        if (n <= CUTOFF)
            return multiplyStandard(A, B);

        // if (n == 1)
        // {
        //     Matrix C(1, std::vector<double>(1));
        //     C[0][0] = A[0][0] * B[0][0];
        //     return C;
        // }

        int mid = n / 2;

        Matrix A11(mid, std::vector<double>(mid));
        Matrix A12(mid, std::vector<double>(mid));
        Matrix A21(mid, std::vector<double>(mid));
        Matrix A22(mid, std::vector<double>(mid));

        Matrix B11(mid, std::vector<double>(mid));
        Matrix B12(mid, std::vector<double>(mid));
        Matrix B21(mid, std::vector<double>(mid));
        Matrix B22(mid, std::vector<double>(mid));

        /* Split matrices */
        for (int i = 0; i < mid; i++)
        {
            for (int j = 0; j < mid; j++)
            {
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j + mid];
                A21[i][j] = A[i + mid][j];
                A22[i][j] = A[i + mid][j + mid];

                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j + mid];
                B21[i][j] = B[i + mid][j];
                B22[i][j] = B[i + mid][j + mid];
            }
        }

        /* 8 recursive multiplications */
        Matrix C11 = add(
            multiplyDivideConquer(A11, B11, mid),
            multiplyDivideConquer(A12, B21, mid));

        Matrix C12 = add(
            multiplyDivideConquer(A11, B12, mid),
            multiplyDivideConquer(A12, B22, mid));

        Matrix C21 = add(
            multiplyDivideConquer(A21, B11, mid),
            multiplyDivideConquer(A22, B21, mid));

        Matrix C22 = add(
            multiplyDivideConquer(A21, B12, mid),
            multiplyDivideConquer(A22, B22, mid));

        /* Combine result */
        Matrix C(n, std::vector<double>(n));
        for (int i = 0; i < mid; i++)
        {
            for (int j = 0; j < mid; j++)
            {
                C[i][j] = C11[i][j];
                C[i][j + mid] = C12[i][j];
                C[i + mid][j] = C21[i][j];
                C[i + mid][j + mid] = C22[i][j];
            }
        }

        return C;
    }
};

#endif