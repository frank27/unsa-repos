#include "DiagProver.h"
#include <fstream>

DiagProver::DiagProver(const std::string &filepath)
{
    std::ifstream f(filepath);
    if (!f)
        throw std::runtime_error("Cannot open: " + filepath);

    f >> n;
    A.resize(n * n);
    for (auto &x : A)
        f >> x;
}

bool DiagProver::prove()
{
    std::cout << "\n=== DIAGONALIZATION THEOREM PROVER ===\n\n";
    printMatrix("A", A);

    eigenDecomposition(); // fills P, eigvals, eigvals_im

    // abort early if complex eigenvalues
    for (int i = 0; i < n; ++i)
        if (std::abs(eigvals_im[i]) > 1e-8)
        {
            printEigenValues();
            std::cout << "\n>>> NOT DIAGONALIZABLE over R (complex eigenvalues)\n\n";
            return false;
        }

    buildD();        // fills D
    buildPInverse(); // fills PInv

    printEigenValues();
    printMatrix("\nP", P);
    printMatrix("\nD", D);
    printMatrix("\nP_inverted", PInv);

    // reconstruct A_reconstructed = P·D·P⁻¹
    auto A_reconstructed = matrixMultiplication(matrixMultiplication(P, D, n), PInv, n);
    printMatrix("\nP*D*P_inverted", A_reconstructed);

    double err = frobeniusError();
    std::cout << "\nA - P*D*P_inverted ^ F = "
              << std::scientific << std::setprecision(4) << err << "\n";

    bool ok = (err < 1e-7);
    std::cout << (ok ? "\nDIAGONALIZABLE  A = P*D*P_inverted  [PROVED]\n\n"
                     : "\nNOT DIAGONALIZABLE  (error exceeds tolerance)\n\n");
    return ok;
}

void DiagProver::eigenDecomposition()
{
    eigvals.resize(n);
    eigvals_im.resize(n, 0.0);
    P.resize(n * n);

    if (isSymmetric())
    {
        // compute eigenvectors of a real symmetric matrix
        // symmetric → dsyev (real eigenvalues guaranteed, P orthogonal)
        std::cout << "[symmetric] using LAPACKE_dsyev\n";
        std::vector<double> tmp = A;
        if (LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n,
                          tmp.data(), n, eigvals.data()) != 0)
            throw std::runtime_error("dsyev failed");
        P = tmp; // ROW_MAJOR output: A[i*n+j] = component i of eigenvector j
    }
    else
    {
        // general → dgeev (may produce complex eigenvalues)
        std::vector<double> tmp = A;
        if (LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n,
                          tmp.data(), n,
                          eigvals.data(), eigvals_im.data(),
                          nullptr, n, P.data(), n) != 0)
            throw std::runtime_error("dgeev failed");
    }
}

void DiagProver::buildD()
{
    D.assign(n * n, 0.0);
    for (int i = 0; i < n; ++i)
        D[i * n + i] = eigvals[i];
}

void DiagProver::buildPInverse()
{
    PInv.resize(n * n);
    if (isSymmetric())
    {
        // orthogonal matrix → P⁻¹ = Pᵀ (free, no LU needed)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                PInv[i * n + j] = P[j * n + i];
    }
    else
    {
        PInv = luInvert(P, n);
    }
}

// ── Frobenius norm ||A - P·D·P⁻¹||_F ────────────────────────────────────────
double DiagProver::frobeniusError() const
{
    auto A_reconstructed = matrixMultiplication(matrixMultiplication(P, D, n), PInv, n);
    double s = 0.0;
    for (int k = 0; k < n * n; ++k)
    {
        double d = A[k] - A_reconstructed[k];
        s += d * d;
    }
    return std::sqrt(s);
}

void DiagProver::printMatrix(const std::string &label,
                             const std::vector<double> &M) const
{
    std::cout << label << " [" << n << "x" << n << "]:\n";
    for (int i = 0; i < n; ++i)
    {
        std::cout << "  ";
        for (int j = 0; j < n; ++j)
        {
            double v = M[i * n + j];
            std::cout << std::setw(11) << std::fixed << std::setprecision(5)
                      << (std::abs(v) < 1e-9 ? 0.0 : v);
        }
        std::cout << "\n";
    }
}

void DiagProver::printEigenValues() const
{
    std::cout << "\nEigenvalues: ";
    for (int i = 0; i < n; ++i)
    {
        std::cout << std::fixed << std::setprecision(5) << eigvals[i];
        if (std::abs(eigvals_im[i]) > 1e-8)
            std::cout << "+" << eigvals_im[i] << "i";
        std::cout << (i < n - 1 ? "  " : "\n");
    }
}

bool DiagProver::isSymmetric() const
{
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (std::abs(A[i * n + j] - A[j * n + i]) > 1e-10)
                return false;
    return true;
}

// static: matrix multiply C = A·B (row-major, n×n)
std::vector<double> DiagProver::matrixMultiplication(const std::vector<double> &A,
                                                     const std::vector<double> &B, int n)
{
    std::vector<double> C(n * n, 0.0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                n, n, n, 1.0, A.data(), n, B.data(), n, 0.0, C.data(), n);
    return C;
}

// invert via LU (dgetrf + dgetri)
std::vector<double> DiagProver::luInvert(std::vector<double> M, int n)
{
    std::vector<int> piv(n);
    if (LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, M.data(), n, piv.data()) != 0)
        throw std::runtime_error("P is singular");
    if (LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, M.data(), n, piv.data()) != 0)
        throw std::runtime_error("dgetri failed");
    return M;
}