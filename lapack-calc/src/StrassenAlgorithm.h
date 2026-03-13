#ifndef STRASSENALGORITHM_H
#define STRASSENALGORITHM_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "Matrix.h"

class StrassenAlgorithm
{
    static constexpr int CUTOFF = 128;

public:
    Matrix multiply(const Matrix &A, const Matrix &B, int n)
    {
        return strassen(A, B, n);
        // Utils::print(a, n, "Strassen Algorithm");
        // return a;
    }

private:
    Matrix add(const Matrix &A, const Matrix &B)
    {
        int n = A.size();
        Matrix C(n, std::vector<double>(n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                C[i][j] = A[i][j] + B[i][j];
        return C;
    }

    Matrix sub(const Matrix &A, const Matrix &B)
    {
        int n = A.size();
        Matrix C(n, std::vector<double>(n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                C[i][j] = A[i][j] - B[i][j];
        return C;
    }

    Matrix multiplyStandard(const Matrix &A, const Matrix &B)
    {
        int n = A.size();
        Matrix C(n, std::vector<double>(n, 0.0));

        for (int i = 0; i < n; ++i)
            for (int k = 0; k < n; ++k)
                for (int j = 0; j < n; ++j)
                    C[i][j] += A[i][k] * B[k][j];

        return C;
    }

    Matrix strassen(const Matrix &A, const Matrix &B, int n)
    {
        // if (n == 1)
        // {
        //     Matrix C(1, std::vector<double>(1));
        //     C[0][0] = A[0][0] * B[0][0];
        //     return C;
        // }
        if (n <= CUTOFF)
        {
            return multiplyStandard(A, B);
        }

        int mid = n / 2;

        Matrix A11(mid, std::vector<double>(mid));
        Matrix A12(mid, std::vector<double>(mid));
        Matrix A21(mid, std::vector<double>(mid));
        Matrix A22(mid, std::vector<double>(mid));

        Matrix B11(mid, std::vector<double>(mid));
        Matrix B12(mid, std::vector<double>(mid));
        Matrix B21(mid, std::vector<double>(mid));
        Matrix B22(mid, std::vector<double>(mid));

        for (int i = 0; i < mid; ++i)
            for (int j = 0; j < mid; ++j)
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

        Matrix M1 = strassen(add(A11, A22), add(B11, B22), mid);
        Matrix M2 = strassen(add(A21, A22), B11, mid);
        Matrix M3 = strassen(A11, sub(B12, B22), mid);
        Matrix M4 = strassen(A22, sub(B21, B11), mid);
        Matrix M5 = strassen(add(A11, A12), B22, mid);
        Matrix M6 = strassen(sub(A21, A11), add(B11, B12), mid);
        Matrix M7 = strassen(sub(A12, A22), add(B21, B22), mid);

        Matrix C11 = add(sub(add(M1, M4), M5), M7);
        Matrix C12 = add(M3, M5);
        Matrix C21 = add(M2, M4);
        Matrix C22 = add(sub(add(M1, M3), M2), M6);

        Matrix C(n, std::vector<double>(n));
        for (int i = 0; i < mid; ++i)
            for (int j = 0; j < mid; ++j)
            {
                C[i][j] = C11[i][j];
                C[i][j + mid] = C12[i][j];
                C[i + mid][j] = C21[i][j];
                C[i + mid][j + mid] = C22[i][j];
            }

        return C;
    }
};

#endif