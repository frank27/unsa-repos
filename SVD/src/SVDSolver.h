#pragma once
//  Decomposes  A (m×n)  →  U (m×m) · Σ (k) · Vᵀ (n×n)
//  where k = min(m,n)
//  Eigenvalue derivation:
//    • eigenvalues  of AᵀA  = σᵢ²   (rows of Vᵀ are the eigenvectors)
//    • eigenvalues  of AAᵀ  = σᵢ²   (columns of U are the eigenvectors)

#include "Matrix.h"
#include <lapacke.h>
#include <cmath>
#include <sstream>

struct SVDResult
{
    Matrix U;                            ///< left  singular vectors (m × m)
    std::vector<double> singular_values; ///< σ₀ ≥ σ₁ ≥ … ≥ 0   (length k)
    Matrix VT;                           ///< Vᵀ: right sing. vectors (n × n)

    std::size_t rank() const noexcept { return singular_values.size(); }

    //  Eigenvalues of AᵀA  (= σᵢ²)
    std::vector<double> eigenvalues_AtA() const
    {
        std::vector<double> ev(singular_values.size());
        for (std::size_t i = 0; i < ev.size(); ++i)
            ev[i] = singular_values[i] * singular_values[i];
        return ev;
    }

    //  Eigenvector i of AᵀA  = row i of Vᵀ
    std::vector<double> eigenvector_AtA(std::size_t i) const
    {
        if (i >= rank())
            throw MatrixError("eigenvector_AtA: index out of range");
        std::size_t n = VT.cols();
        std::vector<double> v(n);
        for (std::size_t j = 0; j < n; ++j)
            v[j] = VT(i, j);
        return v;
    }

    //  Eigenvector i of AAᵀ  = column i of U
    std::vector<double> eigenvector_AAt(std::size_t i) const
    {
        if (i >= rank())
            throw MatrixError("eigenvector_AAt: index out of range");
        std::size_t m = U.rows();
        std::vector<double> u(m);
        for (std::size_t j = 0; j < m; ++j)
            u[j] = U(j, i);
        return u;
    }

    //  Reconstruct  A ≈ U · Σ · Vᵀ  using the first `rank` components
    Matrix reconstruct(std::size_t use_rank = 0) const
    {
        std::size_t m = U.rows();
        std::size_t n = VT.cols();
        std::size_t k = (use_rank == 0) ? rank()
                                        : std::min(use_rank, rank());
        Matrix A(m, n);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t r = 0; r < k; ++r)
                    A(i, j) += U(i, r) * singular_values[r] * VT(r, j);
        return A;
    }

    void print(std::ostream &os = std::cout) const
    {
        os << "\n  Singular values (σ):  [ ";
        for (double s : singular_values)
            os << std::fixed << std::setprecision(6)
               << std::setw(12) << s << " ";
        os << "]\n";
        U.print(os, "U  (left singular vectors)");
        VT.print(os, "Vᵀ (right singular vectors transposed)");

        os << "\n  Eigenvalues of AᵀA  (λᵢ = σᵢ²):\n";
        auto ev = eigenvalues_AtA();
        for (std::size_t i = 0; i < ev.size(); ++i)
            os << "    λ[" << i << "] = "
               << std::fixed << std::setprecision(8) << ev[i] << "\n";

        os << "\n  Eigenvectors of AᵀA (rows of Vᵀ):\n";
        for (std::size_t i = 0; i < rank(); ++i)
        {
            auto v = eigenvector_AtA(i);
            os << "    v[" << i << "] = [ ";
            for (double x : v)
                os << std::fixed << std::setprecision(6)
                   << std::setw(11) << x << " ";
            os << "]\n";
        }

        os << "\n  Eigenvectors of AAᵀ (columns of U):\n";
        for (std::size_t i = 0; i < rank(); ++i)
        {
            auto u = eigenvector_AAt(i);
            os << "    u[" << i << "] = [ ";
            for (double x : u)
                os << std::fixed << std::setprecision(6)
                   << std::setw(11) << x << " ";
            os << "]\n";
        }
    }
};

class SVDSolver
{
public:
    //  Performs A = U Σ Vᵀ via LAPACK dgesdd.
    //    'A' → full U (m×m) and Vᵀ (n×n)   [default]
    static SVDResult compute(const Matrix &A, char jobz = 'A')
    {
        if (A.empty())
            throw MatrixError("SVDSolver::compute — input matrix is empty");
        if (jobz != 'A')
            throw MatrixError("SVDSolver::compute — jobz must be 'A'");

        lapack_int m = static_cast<lapack_int>(A.rows());
        lapack_int n = static_cast<lapack_int>(A.cols());
        lapack_int k = std::min(m, n);

        lapack_int u_cols = (jobz == 'A') ? m : k;
        lapack_int vt_rows = (jobz == 'A') ? n : k;

        Matrix work = A.clone();

        SVDResult result;
        result.singular_values.resize(static_cast<std::size_t>(k));
        result.U = Matrix(static_cast<std::size_t>(m),
                          static_cast<std::size_t>(u_cols));
        result.VT = Matrix(static_cast<std::size_t>(vt_rows),
                           static_cast<std::size_t>(n));

        lapack_int info = LAPACKE_dgesdd(
            LAPACK_ROW_MAJOR,
            jobz,
            m, n,
            work.data(), // overwritten on exit
            n,           // lda  (row-major → lda = cols)
            result.singular_values.data(),
            result.U.data(),
            u_cols, // ldu
            result.VT.data(),
            n // ldvt
        );

        if (info < 0)
        {
            std::ostringstream ss;
            ss << "LAPACKE_dgesdd: illegal argument " << -info;
            throw MatrixError(ss.str());
        }
        if (info > 0)
        {
            throw MatrixError("LAPACKE_dgesdd: DBDSDC did not converge — "
                              "the update process failed (info=" +
                              std::to_string(info) + ")");
        }

        return result;
    }

    //  Frobenius residual  ||A - U Σ Vᵀ||_F
    static double reconstruction_error(const Matrix &A, const SVDResult &res)
    {
        Matrix A_rec = res.reconstruct();
        Matrix diff(A.rows(), A.cols());
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < A.cols(); ++j)
                diff(i, j) = A(i, j) - A_rec(i, j);
        return diff.frobenius_norm();
    }
};