#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <functional>
#include <variant>
#include <string>
#include <iomanip>
#include <iostream>

using Matrix = std::vector<std::vector<double>>;
using ColumnMatrix = std::vector<double>;

struct AnalysisElement
{
    size_t size; // padded power-of-two side length
    Matrix A;    // (Standard / D&C / Strassen)
    Matrix B;
    ColumnMatrix colA; // column-major flat (SGEMM)
    ColumnMatrix colB;
};

using MultiplyAlgorithmDelegate =
    std::variant<
        std::function<void(const Matrix &, const Matrix &, int)>,
        std::function<void(const ColumnMatrix &, const ColumnMatrix &, int)>>;

struct Utils
{
    // Convert row-major Matrix → column-major flat vector.
    // SGEMM convention: element (i,j) lives at index j*n + i.
    static ColumnMatrix toColumnMajor(const Matrix &M, int n)
    {
        ColumnMatrix flat(static_cast<size_t>(n * n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                flat[j * n + i] = M[i][j];
        return flat;
    }

    static void print(const Matrix &A, int n, const std::string &label)
    {
        std::cout << label << "\n"
                  << std::fixed << std::setprecision(15);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
                std::cout << std::setw(30) << A[i][j];
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    static void print(const ColumnMatrix &A, int n, const std::string &label)
    {
        std::cout << label << "\n"
                  << std::fixed << std::setprecision(15);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
                std::cout << std::setw(30) << A[j * n + i];
            std::cout << "\n";
        }
        std::cout << "\n";
    }
};

#endif