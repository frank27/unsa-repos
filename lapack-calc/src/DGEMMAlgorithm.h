#ifndef DGEMMALGORITHM_H
#define DGEMMALGORITHM_H

#include "Matrix.h"
#include <stdexcept>
#include <string>

// ─────────────────────────────────────────────────────────────────────────────
// External DGEMM declaration  (LAPACK/BLAS, Fortran linkage)
// Computes: C = alpha * op(A) * op(B) + beta * C
// ─────────────────────────────────────────────────────────────────────────────
extern "C"
{
    void dgemm_(char *transa, char *transb,
                int *m, int *n, int *k,
                double *alpha,
                double *a, int *lda,
                double *b, int *ldb,
                double *beta,
                double *c, int *ldc);
}

class DGEMMAlgorithm
{
public:
    /**
     * Multiplies two square matrices A and B using LAPACK's DGEMM.
     *
     * DGEMM computes: C = alpha * op(A) * op(B) + beta * C
     * With alpha = 1.0 and beta = 0.0 this simplifies to: C = A * B
     *
     * @param A  Left  input matrix (n×n, double, column-major flat vector)
     * @param B  Right input matrix (n×n, double, column-major flat vector)
     * @param n  Side length of the square matrices
     * @return   Result matrix C = A * B (n×n, double, column-major flat vector)
     */
    static ColumnMatrix Multiply(const ColumnMatrix &A, const ColumnMatrix &B, int n)
    {
        // Output buffer — zero-initialised, column-major, same layout as A and B
        ColumnMatrix C(static_cast<size_t>(n * n), 0.0);

        // ── DGEMM parameters ─────────────────────────────────────────────────
        char transa = 'N';  // Do not transpose A
        char transb = 'N';  // Do not transpose B
        int m = n;          // Rows    of op(A) and C
        int nVal = n;       // Columns of op(B) and C
        int k = n;          // Columns of op(A) = rows of op(B)
        double alpha = 1.0; // Scale factor for A*B  →  C = 1·A·B + 0·C
        int lda = n;        // Leading dimension of A (column-major: lda = rows = n)
        int ldb = n;        // Leading dimension of B
        double beta = 0.0;  // Scale factor for C  (0 = overwrite, not accumulate)
        int ldc = n;        // Leading dimension of C

        // DGEMM has a Fortran-style signature requiring double* not const double*.
        double *pA = const_cast<double *>(A.data());
        double *pB = const_cast<double *>(B.data());
        double *pC = C.data();

        dgemm_(&transa, &transb,
               &m, &nVal, &k,
               &alpha,
               pA, &lda,
               pB, &ldb,
               &beta,
               pC, &ldc);

        return C;
        // Utils::print(C, n, "SGEMM Algorithm");
        // return C;
    }
};

#endif