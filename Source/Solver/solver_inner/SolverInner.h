#ifndef INMOST_SOLVERINNER_H
#define INMOST_SOLVERINNER_H

#include "Source/Solver/SolverInterface.h"
#include "Source/Misc/utils.h"
#include "solver_prototypes.hpp"
#include "solver_bcgsl.hpp"
#define KSOLVER BCGSL_solver

namespace INMOST {

    class SolverInner : public SolverInterface {
    protected:
        Sparse::Matrix *matrix;
        KSOLVER *solver;
        Solver::OrderInfo info;

        INMOST_DATA_ENUM_TYPE maximum_iterations;
        INMOST_DATA_REAL_TYPE atol, rtol, dtol;
    public:
        SolverInner();

        virtual SolverInterface *Copy(const SolverInterface *other);

        virtual void Assign(const SolverInterface *other);

        virtual void Setup(int *argc, char ***argv, SolverParameters &p);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) = 0;

        virtual bool Solve(Sparse::Vector &RHS, Sparse::Vector &SOL);

        virtual bool Clear();

        virtual bool isMatrixSet();

        virtual std::string GetParameter(std::string name) const;

        virtual void SetParameter(std::string name, std::string value);

        virtual INMOST_DATA_ENUM_TYPE Iterations() const;

        virtual INMOST_DATA_REAL_TYPE Residual() const;

        virtual const std::string ReturnReason() const;

        virtual const std::string SolverName() const = 0;

        virtual INMOST_DATA_REAL_TYPE Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxiter);

        virtual void Finalize();

        virtual ~SolverInner();

    };

}


#endif //INMOST_SOLVERINNER_H
