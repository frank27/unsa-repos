#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <map>
#include <functional>
#include <variant>
#include <type_traits>

using Matrix = std::vector<std::vector<double>>;
using ColumnMatrix = std::vector<double>;

struct Utils
{

    static ColumnMatrix toColumnMajor(const Matrix &M, int n)
    {
        ColumnMatrix flat(static_cast<size_t>(n * n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                flat[j * n + i] = M[i][j];
        return flat;
    }

    static Matrix generateMatrix(int size, double minVal, double maxVal, unsigned int seed)
    {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dist(minVal, maxVal);

        Matrix m(size, std::vector<double>(size, 0.0));

        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                m[i][j] = dist(gen);

        return m;
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

    static void WriteToFile(const ColumnMatrix &A, int n, const std::string &filename)
    {
        std::ofstream file(filename);
        if (!file.is_open())
            throw std::runtime_error("Cannot create file: " + filename);

        file << n << "\n";
        file << std::fixed << std::setprecision(15);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                file << std::setw(20) << A[j * n + i];
                if (j < n - 1)
                    file << " ";
            }
            file << "\n";
        }
        file.close();
    }

    static void WriteToFile(const Matrix &A, int n, const std::string &filename)
    {
        std::ofstream file(filename);
        if (!file.is_open())
            throw std::runtime_error("Cannot create file: " + filename);

        file << n << "\n";
        file << std::fixed << std::setprecision(15);

        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                file << std::setw(20) << A[i][j];
                if (j < n - 1)
                    file << " ";
            }
            file << "\n";
            std::cout << std::endl;
        }

        file.close();
    }
};
#endif