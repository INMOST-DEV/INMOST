#ifndef INMOST_SOLVERPETSC_H
#define INMOST_SOLVERPETSC_H

#include "inmost.h"
#include "petsc.h"
#include "solver_petsc.h"


namespace INMOST {

    class SolverPETSc : public SolverInterface {
    private:
		//TODO
		//may be should find another way to count petsc solvers and finalize them
		static unsigned int petscSolversCount;
        std::string parametersFile;
    	KSP* ksp;
    	Mat* matrix;

		INMOST_DATA_ENUM_TYPE local_size, global_size;
    public:
        SolverPETSc();
        SolverPETSc(const SolverInterface* other);

		virtual void Assign(const SolverInterface* other);

        virtual void Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix);
		virtual void SetMatrix(Sparse::Matrix & A, bool ModifiedPattern, bool OldPreconditioner);
        virtual bool Solve(Sparse::Vector &RHS, Sparse::Vector &SOL);
        virtual bool Clear();

    	virtual bool isMatrixSet();

		virtual INMOST_DATA_REAL_TYPE GetPropertyReal(std::string property) const;
		virtual INMOST_DATA_ENUM_TYPE GetPropertyEnum(std::string property) const;

		virtual void SetPropertyReal(std::string property, INMOST_DATA_REAL_TYPE value);
        virtual void SetPropertyEnum(std::string property, INMOST_DATA_ENUM_TYPE value);

		virtual const INMOST_DATA_ENUM_TYPE Iterations() const;
		virtual const INMOST_DATA_REAL_TYPE Residual() const;
		virtual const std::string ReturnReason() const;

		virtual const std::string SolverName() const;

		virtual void Finalize();
        virtual ~SolverPETSc();
    };

}


#endif //INMOST_SOLVERPETSC_H
