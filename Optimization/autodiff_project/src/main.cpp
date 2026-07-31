#include <iostream>
#include "autodiff.h"
#include "Utils.h"

int main() {
    int size = 2;
    int seed = 42;
    double minValue = -100.0;
    double maxValue = 100.0;
    
    Matrix W = Utils::generateMatrix(size, minValue, maxValue, seed);
    Utils::print(W, size, "Matrix W:");
    std::vector<double> x = {1.0, -1.0};

    AutoDiff autodiff(W, x);

    autodiff.forward_and_backward();
    autodiff.printResults();

    return 0;
}
