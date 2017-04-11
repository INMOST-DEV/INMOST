#ifndef INMOST_SOLVERPETSC_H
#define INMOST_SOLVERPETSC_H

#include "petsc.h"
#include "solver_petsc.h"
#include "Source/Solver/SolverInterface.h"
#include "Source/Misc/utils.h"


namespace INMOST {

    class SolverPETSc : public SolverInterface {
    private:
        //TODO
        //may be should find another way to count petsc solvers and finalize them
        static unsigned int petscSolversCount;
        std::string parametersFile;
        KSP *ksp;
        Mat *matrix;

        INMOST_DATA_ENUM_TYPE local_size, global_size;

        INMOST_DATA_ENUM_TYPE maximum_iterations, schwartz_overlap;
        INMOST_DATA_REAL_TYPE atol, dtol, rtol, drop_tolerance, fill_level;
    public:
        SolverPETSc();

        virtual SolverInterface *Copy(const SolverInterface *other);

        virtual void Assign(const SolverInterface *other);

        virtual void Setup(int *argc, char ***argv, SolverParameters &p);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual bool Solve(Sparse::Vector &RHS, Sparse::Vector &SOL);

        virtual bool Clear();

        virtual bool isMatrixSet();

        virtual std::string GetParameter(std::string name) const;

        virtual void SetParameter(std::string name, std::string value);

        virtual INMOST_DATA_ENUM_TYPE Iterations() const;

        virtual INMOST_DATA_REAL_TYPE Residual() const;

        virtual const std::string ReturnReason() const;

        virtual const std::string SolverName() const;

        virtual void Finalize();

        virtual ~SolverPETSc();
    };

}


#endif //INMOST_SOLVERPETSC_H
