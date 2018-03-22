#include "SolverInterface.h"

#if defined(USE_SOLVER)
namespace INMOST {

    SolverParameters::SolverParameters(std::string solverName, std::string solverPrefix, std::string internalFile)
            : solverName(solverName),
              solverPrefix(solverPrefix),
              internalFile(internalFile) {}

    SolverParameters::SolverParameters(const SolverParameters &other) : solverName(other.solverName),
                                                                        solverPrefix(other.solverPrefix),
                                                                        internalFile(other.internalFile),
                                                                        parameters(other.parameters) {
    }


    SolverParameters::~SolverParameters() {}

    void SolverParameters::SetInnerParametersFromFile(const std::string &file, SolverInterface *solver) {
        char line[4096];
        char parameterName[4096];
        char parameterValue[4096];
        FILE *databaseFile = fopen(file.c_str(), "r");
        if (!databaseFile) return;
        while (!feof(databaseFile) && fgets(line, 4096, databaseFile)) {
            if (line[0] == '#') continue;
            sscanf(line, "%s %s", parameterName, parameterValue);
            solver->SetParameter(parameterName, parameterValue);
        }
    }

}
#endif //USE_SOLVER