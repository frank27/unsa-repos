#pragma once

#include <stdexcept>
#include <string>
#include "Utils.h"
#include <lapacke.h>
#include <cblas.h>

class DGEMMRoutine
{
public:
    /**
     * Multiplies two square matrices A and B using CBLAS's cblas_dgemm.
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
        // Output buffer — zero-initialised, column-major
        ColumnMatrix C(static_cast<size_t>(n * n), 0.0);

        cblas_dgemm(
            CblasColMajor,  // Column-major storage order
            CblasNoTrans,   // Do not transpose A
            CblasNoTrans,   // Do not transpose B
            n,              // Rows    of op(A) and C
            n,              // Columns of op(B) and C
            n,              // Columns of op(A) = rows of op(B)
            1.0,            // alpha: scale factor for A*B  →  C = 1·A·B + 0·C
            A.data(),       // Pointer to A (const-correct, no cast needed)
            n,              // Leading dimension of A (column-major: lda = rows = n)
            B.data(),       // Pointer to B
            n,              // Leading dimension of B
            0.0,            // beta: scale factor for C (0 = overwrite, not accumulate)
            C.data(),       // Pointer to C (output)
            n               // Leading dimension of C
        );

        return C;
    }
};