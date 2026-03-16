#include <iostream>
#include "DGEMMRoutine.h"
#include "StandardAlgorithm.h"
#include "Utils.h"

void clearFiles(std::string file1, std::string file2, std::string fileOutSGEMM, std::string fileOutStandard)
{
    std::ofstream(file1, std::ios::trunc).close();
    std::ofstream(file2, std::ios::trunc).close();
    std::ofstream(fileOutSGEMM, std::ios::trunc).close();
    std::ofstream(fileOutStandard, std::ios::trunc).close();
}

int main()
{
    try
    {
        // DATA_DIR is defined by CMake and points to the src/ folder.
        const std::string dataDir = DATA_DIR;

        std::string file1 = dataDir + "matrix_a.txt";
        std::string file2 = dataDir + "matrix_b.txt";
        std::string fileOutSGEMM = dataDir + "resultDGEMM.txt";
        std::string fileOutStandard = dataDir + "resultStandard.txt";
        clearFiles(file1, file2, fileOutSGEMM, fileOutStandard);

        const unsigned int SEED_A = 4;
        const unsigned int SEED_B= 10;
        const double MIN_VAL = -10.0;
        const double MAX_VAL = 10.0;
        const unsigned int SIZE = 4;

        Matrix A_2D = Utils::generateMatrix(SIZE, MIN_VAL, MAX_VAL, SEED_A);
        Matrix B_2D = Utils::generateMatrix(SIZE, MIN_VAL, MAX_VAL, SEED_B);

        Utils::WriteToFile(A_2D, SIZE, file1);
        Utils::WriteToFile(B_2D, SIZE, file2);

        ColumnMatrix A = Utils::toColumnMajor(A_2D, SIZE);
        ColumnMatrix B = Utils::toColumnMajor(B_2D, SIZE);

        ColumnMatrix C = DGEMMRoutine::Multiply(A, B, SIZE);
        Utils::print(C, SIZE, "SGEMM algorithm: Matrix C = A x B (result):");
        Utils::WriteToFile(C, SIZE, fileOutSGEMM);

        auto C2D = StandardAlgorithm::Multiply(A_2D, B_2D, SIZE);
        Utils::print(C2D, SIZE, "Standard algorithm: Matrix C = A x B (result):");
        Utils::WriteToFile(C2D, SIZE, fileOutStandard);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

