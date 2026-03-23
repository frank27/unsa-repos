#ifndef LINEAR_ALGEBRA_ASSIGNMENTS_H
#define LINEAR_ALGEBRA_ASSIGNMENTS_H

#include "Utils.h"
#include <vector>
#include <lapacke.h>

class LinearAlgebraAssignments
{
public:
    void execute()
    {
        assignments1();
        assignments2();
        assignments3();
        assignments4();
    }

private:
    // 1) Cálculo de la Norma Euclidiana - norma de Frobenius
    double calculateEuclideanNorm(const std::vector<double> &vector)
    {
        char norm = 'F';
        int m = 1;
        int n = vector.size();
        const double * a = vector.data();
        int32_t lda = 1;

        return LAPACKE_dlange(LAPACK_COL_MAJOR, norm, n, m, a, lda);
    }

    void assignments1()
    {
        std::vector<double> vector = {3.0, 4.0};
        double norm = calculateEuclideanNorm(vector);
        std::cout << "Tarea 1:  La norma euclidiana del vector es: " << norm << std::endl;
        std::cout << std::endl;
    }

    // 2) Determinante de una Matriz
    double calculateDeterminant(std::vector<double> &matrix, int n)
    {
        int lda = n;
        std::vector<int> ipiv(n);
        int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, matrix.data(), lda, ipiv.data());

        if (info != 0)
        {
            return 0.0; // singular matrix
        }

        double det = 1.0;
        int swaps = 0;

        for (int i = 0; i < n; i++)
        {
            det *= matrix[i * n + i]; // product of diagonal elements

            if (ipiv[i] != i + 1)
            {
                swaps++;
            }
        }

        if (swaps % 2 != 0)
            det *= -1.0;

        return det;
    }

    void assignments2()
    {
        std::vector<double> matrix = {
            1, 2, 3,
            4, 5, 6,
            7, 8, 9};
        int n = 3;

        double det = calculateDeterminant(matrix, n);
        std::cout << "Tarea 2: La determinante de la matriz es: " << det << std::endl;
        std::cout << std::endl;
    }

    // 3) Determinación del Rango
    int calculateRank(std::vector<double> &matrix, int m, int n)
    {
        int k = std::min(m, n);
        std::vector<double> S(k);
        std::vector<double> U(m * m);
        std::vector<double> VT(n * n);
        std::vector<double> superb(k - 1);
        int lda = n;

        int info = LAPACKE_dgesvd(
            LAPACK_ROW_MAJOR,
            'A',
            'A',
            m, n,
            matrix.data(),
            lda,
            S.data(),
            U.data(),
            m,
            VT.data(),
            n,
            superb.data());

        if (info != 0)
            return -1;

        int rank = 0;
        for (double s : S)
        {
            if (s > 1e-10)
                rank++;
        }
        return rank;
    }

    void assignments3()
    {
        std::vector<double> matrix = {
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0};
        int m = 3, n = 3;

        int rank = calculateRank(matrix, m, n);
        std::cout << "Tarea 3: El rango de la matriz es: " << rank << std::endl;
        std::cout << std::endl;
    }

    // 4) Cálculo de la Inversa de una Matriz
    std::tuple<bool, std::string> invertMatrix(std::vector<double> &matrix, int n)
    {
        std::vector<int> ipiv(n);
        int lda = n, info;

        // LU factorization
        info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, matrix.data(), lda, ipiv.data());
        if (info != 0)
        {
            if (info > 0)
                return std::make_tuple(false, "Error: matrix is singular at U(" + std::to_string(info) + "," + std::to_string(info) + "), cannot invert.");

            return std::make_tuple(false, "Error: invalid argument " + std::to_string(-info) + " in dgetrf.");
        }

        info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, matrix.data(), lda, ipiv.data());
        if (info != 0)
        {
            if (info > 0)
                return std::make_tuple(false, "Error: matrix is singular at U(" + std::to_string(info) + "," + std::to_string(info) + "), cannot invert.");
            return std::make_tuple(false, "Error: invalid argument " + std::to_string(-info) + " in dgetri.");
        }

        return std::make_tuple(true, "");
    }

    void assignments4()
    {
        std::vector<double> matrix = {
            3.0, -1.0, -1.0,
            2.0, 1.0, 0.0,
            3.0, 1.0, 2.0};
        int n = 3;

        std::tuple<bool, std::string> result = invertMatrix(matrix, n);

        if (std::get<0>(result))
        {
            Utils::print(matrix, n, "Tarea 4: La inverza de matriz es:");
        }
        else
        {
            std::cout << std::get<1>(result) << std::endl;
        }
    }
};

#endif