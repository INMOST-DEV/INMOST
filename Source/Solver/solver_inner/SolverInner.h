#ifndef INMOST_SOLVERINNER_H
#define INMOST_SOLVERINNER_H

#include "inmost_solver_interface.h"
#include "solver_prototypes.hpp"
#include "solver_bcgsl.hpp"


namespace INMOST {

    class SolverInner : public SolverInterface {
    protected:
        Sparse::Matrix *matrix;
        KSOLVER *solver;
        Solver::OrderInfo info;
    public:
        SolverInner(SolverParameters &parameters);

        SolverInner(const SolverInterface *other);

        virtual void Assign(const SolverInterface *other);

        virtual void Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) = 0;

        virtual bool Solve(Sparse::Vector &RHS, Sparse::Vector &SOL);

        virtual bool Clear();

        virtual bool isMatrixSet();

        virtual const INMOST_DATA_ENUM_TYPE Iterations() const;

        virtual const INMOST_DATA_REAL_TYPE Residual() const;

        virtual const std::string ReturnReason() const;

        virtual const std::string SolverName() const = 0;

        virtual const INMOST_DATA_REAL_TYPE Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxiter);

        virtual void Finalize();

        virtual ~SolverInner();

    };

}


#endif //INMOST_SOLVERINNER_H
