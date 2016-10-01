#include "SolverInner.h"

namespace INMOST {

    SolverInner::SolverInner() {

    }

    SolverInner::SolverInner(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverInner::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverInner::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {
        FILE *databaseFile = fopen(parameters_file, "r");
        if (!databaseFile) {
            std::cout << "Inner options file not found" << std::endl;
            return;
        }
        char *tmp = (char *) calloc(256, sizeof(char));
        char *parameterName = (char *) calloc(128, sizeof(char));
        char *parameterValue = (char *) calloc(128, sizeof(char));
        char *type = (char *) calloc(36, sizeof(char));
        while (!feof(databaseFile) && fgets(tmp, 256, databaseFile)) {
            //removeSpaces(tmp);
            char *line = tmp;
            //Comment str
            if (line[0] == '#') continue;
            //First 4 chars is 'real' or 'enum'
            bool isReal = false, isEnum = false;
            sscanf(line, "%s %s %s", type, parameterName, parameterValue);
            if (strncmp(type, "real", 4) == 0) {
                //fprintf(stdout, "Real parameter:");
                isReal = true;
            } else if (strncmp(type, "enum", 4) == 0) {
                //fprintf(stdout, "Enum parameter:");
                isEnum = true;
            } else {
                fprintf(stderr, "Skipping bad line: %s", line);
                continue;
            }
            if (isReal) {
                parameters.SetParameterReal(parameterName, atof(parameterValue));
            } else {
                parameters.SetParameterEnum(parameterName, static_cast<INMOST_DATA_ENUM_TYPE>(atoi(parameterValue)));
            }
        }
        free(parameterValue);
        free(parameterName);
        free(tmp);
        free(type);
    }

    bool SolverInner::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = parameters.GetParameterEnum("maximum_iterations");
        solver->RealParameter("rtol") = parameters.GetParameterReal("relative_tolerance");
        solver->RealParameter("atol") = parameters.GetParameterReal("absolute_tolerance");
        solver->RealParameter("divtol") = parameters.GetParameterReal("divergence_tolerance");

        return solver->Solve(RHS, SOL);
    }

    bool SolverInner::Clear() {
        info.Clear();
        if (matrix != NULL) {
            delete matrix;
            matrix = NULL;
        }
        if (solver != NULL) {
            delete solver;
            solver = NULL;
        }
        return true;
    }

    bool SolverInner::isMatrixSet() {
        return matrix != NULL;
    }

    const INMOST_DATA_ENUM_TYPE SolverInner::Iterations() const {
        return solver->GetIterations();
    }

    const INMOST_DATA_REAL_TYPE SolverInner::Residual() const {
        return solver->GetResidual();
    }

    const std::string SolverInner::ReturnReason() const {
        return solver->GetReason();
    }

    void SolverInner::Finalize() {

    }

    SolverInner::~SolverInner() {
        this->Clear();
    }
}