#include "Diagonalization.h"
#include "DiagProver.h"

void Diagonalization::Execute()
{
    const std::string dataDir = DATA_DIR;
    try
    {
        const std::string path = dataDir+"DiagData/matrix.txt";
        DiagProver prover(path);
        prover.prove();
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
    }
}