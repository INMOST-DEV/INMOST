#ifndef INMOST_SOLVERHYPRE_H
#define INMOST_SOLVERHYPRE_H

#include "../Solver/SolverInterface.h"

//#include <HYPRE.h>
#if !defined(USE_MPI)
#define HYPRE_SEQUENTIAL
#endif
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

namespace INMOST {

    enum HypreType { BoomerAMG, ParaSails, AMS, PILUT, Euclid};

    class SolverHYPRE : public SolverInterface 
    {
    private:
        HYPRE_ParCSRMatrix parcsr_A;
        HYPRE_Solver solver, precond;
        mutable HYPRE_IJVector b, x;
        mutable HYPRE_ParVector par_b, par_x;
        mutable std::vector<int> rows;
        int size;
        HypreType type;
    public:
        SolverHYPRE(HypreType type);

        virtual SolverInterface *Copy(const SolverInterface *other);

        virtual void Assign(const SolverInterface *other);

        virtual void Setup(int *argc, char ***argv, SolverParameters &p);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

        virtual bool Clear();

        virtual bool isMatrixSet();

        virtual std::string GetParameter(std::string name) const;

        virtual void SetParameter(std::string name, std::string value);

        virtual INMOST_DATA_ENUM_TYPE Iterations() const;

        virtual INMOST_DATA_REAL_TYPE Residual() const;

        virtual const std::string ReturnReason() const;

        virtual const std::string SolverName() const;

        virtual void Finalize();

        virtual ~SolverHYPRE();

    };

}


#endif //INMOST_SOLVERSUPERLU_H
