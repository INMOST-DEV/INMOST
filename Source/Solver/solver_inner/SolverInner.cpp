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
            return;
        }
        char *tmp = (char *) calloc(256, sizeof(char));
        char *parameterName = (char *) calloc(128, sizeof(char));
        char *parameterValue = (char *) calloc(128, sizeof(char));
        while (!feof(databaseFile) && fgets(tmp, 256, databaseFile)) {
            char *line = tmp;
            //Comment str
            if (line[0] == '#') continue;
            sscanf(line, "%s %s", parameterName, parameterValue);
            parameters.SetParameter(parameterName, parameterValue);
        }
        free(parameterValue);
        free(parameterName);
        free(tmp);
    }

    bool SolverInner::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = parameters.get<INMOST_DATA_ENUM_TYPE>("maximum_iterations");
        solver->RealParameter("rtol") = parameters.get<INMOST_DATA_REAL_TYPE>("relative_tolerance");
        solver->RealParameter("atol") = parameters.get<INMOST_DATA_REAL_TYPE>("absolute_tolerance");
        solver->RealParameter("divtol") = parameters.get<INMOST_DATA_REAL_TYPE>("divergence_tolerance");

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

    void SolverInner::SetDefaultParameters() {
        this->SetParameter("additive_schwartz_overlap", "1");
        this->SetParameter("maximum_iterations", "2500");
        this->SetParameter("reorder_nonzero", "1");
        this->SetParameter("rescale_iterations", "6");
        this->SetParameter("condition_estimation", "1");
        this->SetParameter("adapt_ddpq_tolerance", "1");
        this->SetParameter("gmres_substeps", "2");

        this->SetParameter("absolute_tolerance", "1.0e-5");
        this->SetParameter("relative_tolerance", "1.0e-12");
        this->SetParameter("divergence_tolerance", "1.0e+100");
        this->SetParameter("drop_tolerance", "0.005");
        this->SetParameter("reuse_tolerance", "0.00005");
        this->SetParameter("ddpq_tolerance", "0.75");
        this->SetParameter("fill_level", "3");
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