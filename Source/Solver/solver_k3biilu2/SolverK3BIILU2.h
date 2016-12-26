#ifndef INMOST_SOLVERK3BIILU2_H
#define INMOST_SOLVERK3BIILU2_H

#include "Source/Solver/SolverInterface.h"
#include "solver_k3biilu2.h"

namespace INMOST {

    class SolverK3BIILU2 : public SolverInterface {
    private:
        bcg_k3biilu2 *solver_data;
        matrix_k3biilu2 *matrix_data;
        INMOST_DATA_ENUM_TYPE local_size, global_size;
    public:
        SolverK3BIILU2();

        virtual SolverInterface *Copy(const SolverInterface *other);

        virtual void Assign(const SolverInterface *other);

        virtual void Setup(int *argc, char ***argv, SolverParameters &p);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

        virtual bool Clear();

        virtual bool isMatrixSet();

        virtual std::string GetParameter(std::string name) const;

        virtual void SetParameter(std::string name, std::string value);

        virtual const INMOST_DATA_ENUM_TYPE Iterations() const;

        virtual const INMOST_DATA_REAL_TYPE Residual() const;

        virtual const std::string ReturnReason() const;

        virtual const std::string SolverName() const;

        virtual void Finalize();

        virtual ~SolverK3BIILU2();

    };

}


#endif //INMOST_SOLVERK3BIILU2_H
