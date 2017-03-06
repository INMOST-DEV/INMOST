#include "SolverInterface.h"
#if defined(USE_SOLVER)
namespace INMOST {

    SolverParameters::SolverParameters(std::string solverName, std::string solverPrefix, std::string internalFile) : solverName(solverName),
                                                                                                                             solverPrefix(solverPrefix),
                                                                                                                             internalFile(internalFile) {}

    SolverParameters::SolverParameters(const SolverParameters &other) : solverName(other.solverName), solverPrefix(other.solverPrefix),
                                                                                internalFile(other.internalFile),
                                                                                parameters(other.parameters) {
    }


    SolverParameters::~SolverParameters() {}

}
#endif //USE_SOLVER