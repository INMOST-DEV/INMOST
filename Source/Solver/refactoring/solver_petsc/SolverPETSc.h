//
// Created by Dmitri Bagaev on 22.09.16.
//

#ifndef INMOST_SOLVERPETSC_H
#define INMOST_SOLVERPETSC_H

#include "petsc.h"
#include "inmost.h"
#include "new_solver_petsc.h"
#include "Source/Solver/solver_petsc.h"
#include "Source/Solver/refactoring/SolverFactory.h"
#include "Source/Solver/refactoring/SolverInterface.h"


namespace INMOST {

    class SolverPETSc : public SolverInterface {
    private:
        std::string parameters_file;
    	KSP* ksp;
    	Mat* matrix;
        Vec* rhs;
		Vec* solution;

		INMOST_DATA_ENUM_TYPE local_size, global_size;
    public:
        SolverPETSc();
        SolverPETSc(const SolverInterface* other);
        virtual const std::string getSolverName() const;
        virtual void Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix);
        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);
        virtual void SetMatrix(Sparse::Matrix & A, bool ModifiedPattern, bool OldPreconditioner);
        virtual void Assign(const SolverInterface* other);
        virtual const INMOST_DATA_ENUM_TYPE Iterations() const;
    	virtual const INMOST_DATA_REAL_TYPE Residual() const;
    	virtual const std::string ReturnReason() const;
    	virtual bool isMatrixSet();
        virtual void Finalize();
        virtual ~SolverPETSc();
    };

}


#endif //INMOST_SOLVERPETSC_H
