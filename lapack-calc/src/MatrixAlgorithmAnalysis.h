#ifndef MATRIXALGORITHMANALYSIS_H
#define MATRIXALGORITHMANALYSIS_H

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <map>
#include <functional>
#include <variant>
#include <type_traits>

#include "Matrix.h"
#include "StandardAlgorithm.h"
#include "DivideConquerAlgorithm.h"
#include "StrassenAlgorithm.h"
#include "DGEMMAlgorithm.h"

enum class AlgorithmType
{
    Standard = 0,
    DivideAndConquer = 1,
    Strassen = 2,
    DGEMM = 3
};

class MatrixAlgorithmAnalysis
{
public:
    // const std::vector<int> SIZES = {2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
    // const std::vector<int> SIZES = {2, 4, 8, 16, 32, 64, 128, 256, 512};
    const std::vector<int> SIZES = {2, 4, 8, 16};
    // const std::vector<int> SIZES = {4};
    const int REPEATS = 3;

    std::vector<std::tuple<AlgorithmType, int, double, int>> results;

    const std::map<AlgorithmType, std::string> ALGORITHM_LEGEND = {
        {AlgorithmType::Standard, "Standard Algorithm"},
        {AlgorithmType::DivideAndConquer, "Divide And Conquer"},
        {AlgorithmType::Strassen, "Strassen Algorithm"},
        {AlgorithmType::DGEMM, "DGEMM Algorithm"},
    };

    void execute()
    {
        const auto dataset = buildDataset();
        const auto methods = getAllMethods();

        for (const auto &data : dataset)
        {
            for (const auto &[id, delegate] : methods)
            {
                for (int rep = 0; rep < REPEATS; ++rep)
                {
                    const double t = measureTime(delegate, data);
                    results.emplace_back(id, static_cast<int>(data.size), t, rep);
                    std::cout << "Algo " << static_cast<int>(id)
                              << " Size " << data.size
                              << " Time " << t << "s\n";
                }
            }
        }

        exportData();
    }

private:
    // ===================== ALGORITHMS =====================
    std::map<AlgorithmType, MultiplyAlgorithmDelegate> getAllMethods()
    {
        static StandardAlgorithm standard;
        static DivideConquerAlgorithm divideConquer;
        static StrassenAlgorithm strassen;
        static DGEMMAlgorithm sgemm;

        using RowFn = std::function<void(const Matrix &, const Matrix &, int)>;
        using ColFn = std::function<void(const ColumnMatrix &, const ColumnMatrix &, int)>;

        return {
            {AlgorithmType::Standard,
             RowFn{[&](const Matrix &A, const Matrix &B, int n)
                   { standard.multiply(A, B, n); }}},
            {AlgorithmType::Strassen,
             RowFn{[&](const Matrix &A, const Matrix &B, int n)
                   { strassen.multiply(A, B, n); }}},

            {AlgorithmType::DGEMM,
             ColFn{[&](const ColumnMatrix &A, const ColumnMatrix &B, int n)
                   { sgemm.Multiply(A, B, n); }}},
        };
    }

    // ===================== DATASET =====================
    std::vector<AnalysisElement> buildDataset()
    {
        std::vector<AnalysisElement> dataset;
        dataset.reserve(SIZES.size());

        for (const int size : SIZES)
        {
            Matrix rowA = generateRandomMatrix(size, size, 0.0, 100.0);
            Matrix rowB = generateRandomMatrix(size, size, 0.0, 100.0);
            const int n = static_cast<int>(rowA.size()); // padded power-of-two

            dataset.push_back({static_cast<size_t>(n),
                               rowA,
                               rowB,
                               Utils::toColumnMajor(rowA, n),
                               Utils::toColumnMajor(rowB, n)});

            // Utils::print(rowA, n, "Matrix A");
            // Utils::print(rowB, n, "Matrix B");

            // Utils::print(Utils::toColumnMajor(rowA, n), n, "Matrix Column A");
            // Utils::print(Utils::toColumnMajor(rowB, n), n, "Matrix Column B");
        }
        return dataset;
    }

    // ===================== TIMING =====================

    // std::visit dispatches to the correct overload at runtime:
    //   - RowFn  → uses data.A   / data.B   (row-major Matrix)
    //   - ColFn  → uses data.colA / data.colB (column-major ColumnMatrix)
    double measureTime(const MultiplyAlgorithmDelegate &delegate,
                       const AnalysisElement &data)
    {
        const int n = static_cast<int>(data.size);

        const auto start = std::chrono::high_resolution_clock::now();

        std::visit([&](const auto &fn)
                   {
            using Fn = std::decay_t<decltype(fn)>;
            if constexpr (std::is_invocable_v<Fn, const Matrix &, const Matrix &, int>)
                fn(data.A, data.B, n);
            else
                fn(data.colA, data.colB, n); }, delegate);

        const auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start).count();
    }

    // ===================== UTILS =====================

    int nextPowerOfTwo(int n)
    {
        int p = 1;
        while (p < n)
            p <<= 1;
        return p;
    }

    Matrix generateRandomMatrix(int rows, int cols, double minVal, double maxVal)
    {
        std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<double> dist(minVal, maxVal);

        const int size = nextPowerOfTwo(std::max(rows, cols));
        Matrix m(size, std::vector<double>(size, 0.0));

        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                m[i][j] = dist(gen);

        return m;
    }

    void exportData()
    {
        std::ofstream file("result_data.csv");
        file << "Algorithm,Size,RepeatNumber,TimeSeconds,AlgorithmName\n";

        for (const auto &[algo, size, time, rep] : results)
            file << static_cast<int>(algo) << ","
                 << size << ","
                 << rep << ","
                 << time << ","
                 << ALGORITHM_LEGEND.at(algo) << "\n";

        std::cout << "CSV exported successfully.\n";
    }
};

#endif