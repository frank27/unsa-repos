#ifndef AUTO_DIFF_H
#define AUTO_DIFF_H

#include <vector>
#include <iostream>
#include <cblas.h>  // Para operaciones BLAS
#include "Utils.h"

class AutoDiff {
private:
    Matrix W;  // Matriz de coeficientes
    std::vector<double> x;               // Vector de entrada
    std::vector<double> y;               // Resultado de f(Wx)
    std::vector<double> grad;            // Gradiente df/dx

public:
    AutoDiff(const Matrix& W, const std::vector<double>& x);
    void forward_and_backward();
    void printResults();
};

#endif
