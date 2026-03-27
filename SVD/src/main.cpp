#include "Matrix.h"
#include "SVDSolver.h"

#include <iostream>
#include <string>
#include <vector>

static void runSVD(const Matrix &A)
{
    A.print(std::cout, "Input A");

    SVDResult res = SVDSolver::compute(A);
    res.print();

    double err = SVDSolver::reconstruction_error(A, res);
    double tolerance = 1e-15;
    std::cout << "\n  ||A - U·Σ·Vᵀ||_F = "
              << std::scientific << std::setprecision(4) << err << "\n";
    if (err < tolerance)
    {
        std::cout << "Valid reconstruction (||A - U·Σ·Vᵀ||_F < " << tolerance << ")\n";
    }
    else
    {
        std::cout << "Invalid reconstruction (||A - U·Σ·Vᵀ||_F > " << tolerance << ")\n";
    }
}

int main(int argc, char *argv[])
{
    const std::string dataDir = DATA_DIR;

    std::string file_path = dataDir + "dataset/matrix.txt";

    try
    {
        if (!file_path.empty())
        {
            Matrix A = Matrix::from_file(file_path);
            runSVD(A);
        }
    }
    catch (const MatrixError &e)
    {
        std::cerr << "\n"
                  << e.what() << "\n";
        return 1;
    }
    catch (const std::exception &e)
    {
        std::cerr << "\n[std::exception] " << e.what() << "\n";
        return 2;
    }
    return 0;
}