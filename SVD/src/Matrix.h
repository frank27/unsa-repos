#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

class MatrixError : public std::runtime_error
{
public:
    explicit MatrixError(const std::string &msg)
        : std::runtime_error("[MatrixError] " + msg) {}
};

class Matrix
{
public:
    Matrix() : numRows(0), numCols(0) {}

    Matrix(std::size_t rows, std::size_t cols)
        : numRows(rows), numCols(cols), vData(rows * cols, 0.0)
    {
        if (rows == 0 || cols == 0)
            throw MatrixError("Dimensions must be > 0");
    }

    /// Construct from a flat row-major vector
    Matrix(std::size_t rows, std::size_t cols, std::vector<double> data)
        : numRows(rows), numCols(cols), vData(std::move(data))
    {
        if (vData.size() != numRows * numCols)
            throw MatrixError("Data size does not match dimensions");
    }
    double &operator()(std::size_t i, std::size_t j)
    {
        bounds_check(i, j);
        return vData[i * numCols + j];
    }
    double operator()(std::size_t i, std::size_t j) const
    {
        bounds_check(i, j);
        return vData[i * numCols + j];
    }

    double *data() noexcept { return vData.data(); }
    const double *data() const noexcept { return vData.data(); }

    std::size_t rows() const noexcept { return numRows; }
    std::size_t cols() const noexcept { return numCols; }
    std::size_t size() const noexcept { return vData.size(); }
    bool empty() const noexcept { return vData.empty(); }

    static Matrix from_file(const std::string &path)
    {
        std::ifstream f(path);
        if (!f.is_open())
            throw MatrixError("Cannot open file: " + path);

        std::size_t rows = 0, cols = 0;

        // Skip comment lines and read the header
        std::string line;
        while (std::getline(f, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            std::istringstream ss(line);
            if (!(ss >> rows >> cols))
                throw MatrixError("Bad header in " + path +
                                  " — expected: rows cols");
            break;
        }
        if (rows == 0 || cols == 0)
            throw MatrixError("Header not found or zero dimensions in " + path);

        std::vector<double> data;
        data.reserve(rows * cols);

        while (std::getline(f, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            std::istringstream ss(line);
            double v;
            while (ss >> v)
                data.push_back(v);
        }

        if (data.size() != rows * cols)
        {
            std::ostringstream err;
            err << "Expected " << rows * cols
                << " values but found " << data.size()
                << " in " << path;
            throw MatrixError(err.str());
        }

        return Matrix(rows, cols, std::move(data));
    }

    void to_file(const std::string &path, int precision = 8) const
    {
        std::ofstream f(path);
        if (!f.is_open())
            throw MatrixError("Cannot write file: " + path);

        f << "# Matrix written by Matrix::to_file\n";
        f << numRows << " " << numCols << "\n";
        for (std::size_t i = 0; i < numRows; ++i)
        {
            for (std::size_t j = 0; j < numCols; ++j)
            {
                f << std::setw(precision + 6) << std::fixed
                  << std::setprecision(precision)
                  << vData[i * numCols + j];
                if (j + 1 < numCols)
                    f << "  ";
            }
            f << "\n";
        }
    }

    void print(std::ostream &os = std::cout,
               const std::string &name = "",
               int precision = 6) const
    {
        if (!name.empty())
            os << "\n  [" << name << "]  (" << numRows << " x " << numCols << ")\n";
        for (std::size_t i = 0; i < numRows; ++i)
        {
            os << "  | ";
            for (std::size_t j = 0; j < numCols; ++j)
                os << std::setw(precision + 5) << std::fixed
                   << std::setprecision(precision)
                   << vData[i * numCols + j] << " ";
            os << "|\n";
        }
    }

    //  Frobenius norm  ||A||_F
    double frobenius_norm() const noexcept
    {
        double s = 0.0;
        for (double v : vData)
            s += v * v;
        return std::sqrt(s);
    }

    Matrix clone() const { return *this; }

private:
    std::size_t numRows;
    std::size_t numCols;
    std::vector<double> vData;

    void bounds_check(std::size_t i, std::size_t j) const
    {
        if (i >= numRows || j >= numCols)
        {
            std::ostringstream ss;
            ss << "Index (" << i << "," << j << ") out of range for "
               << numRows << "x" << numCols << " matrix";
            throw MatrixError(ss.str());
        }
    }
};