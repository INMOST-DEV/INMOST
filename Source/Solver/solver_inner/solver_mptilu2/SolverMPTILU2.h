#ifndef INMOST_SOLVERMPTILU2_H
#define INMOST_SOLVERMPTILU2_H


#include <inmost.h>
#include "solver_mtilu2.hpp"
#include "../solver_bcgsl.hpp"

namespace INMOST {

    class SolverMPTILU2 : public SolverInterface {
    private:
        Sparse::Matrix *matrix;
        KSOLVER *solver;
        Solver::OrderInfo info;
    public:
        SolverMPTILU2();

        SolverMPTILU2(const SolverInterface *other);

        virtual void Assign(const SolverInterface *other);

        virtual void Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual bool Solve(Sparse::Vector &RHS, Sparse::Vector &SOL);

        virtual bool Clear();

        virtual bool isMatrixSet();

        virtual const INMOST_DATA_ENUM_TYPE Iterations() const;

        virtual const INMOST_DATA_REAL_TYPE Residual() const;

        virtual const std::string ReturnReason() const;

        virtual const std::string SolverName() const;

        virtual void Finalize();

        virtual ~SolverMPTILU2();
    };

}


#endif //INMOST_SOLVERMPTILU2_H
