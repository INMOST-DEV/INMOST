#ifndef INMOST_SOLVERK3BIILU2_H
#define INMOST_SOLVERK3BIILU2_H

#include "inmost_solver_interface.h"
#include "solver_k3biilu2.h"

namespace INMOST {

    class SolverK3BIILU2 : public SolverInterface {
    private:
        void *solver_data;
        void *matrix_data;
        INMOST_DATA_ENUM_TYPE local_size, global_size;
    public:
        SolverK3BIILU2();

        SolverK3BIILU2(const SolverInterface *other);

        virtual void Assign(const SolverInterface *other);

        virtual void Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

        virtual bool Clear();

        virtual bool isMatrixSet();

        virtual void SetDefaultParameters();

        virtual SolverParameter GetParameter(std::string name) const;

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
